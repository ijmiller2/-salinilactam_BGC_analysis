#!/usr/bin/env python

import pandas as pd
import argparse
import subprocess
import os
import re

#1. Use dict to select top blasthit, based on eval
#2. Load pathway as dict of dict
#3. Report alignment info

parser = argparse.ArgumentParser(description='Script to calculate BGC fragmentation\
    based on a blastn alignment')
parser.add_argument('-a','--assembly', help='De novo assembly nucl fasta file', required=True)
parser.add_argument('-r','--reference', help='Reference assembly nucl fasta file', required=True)
parser.add_argument('-p','--pathway_table', help='Tab separated pathway table\
    with min and max nucleotide coordinates', required=True)
args = vars(parser.parse_args())

de_novo_assembly = args['assembly']
blastn_output = os.path.splitext(de_novo_assembly)[0] + ".blastn"
reference_assembly = args['reference']

def run_blastn(query_path,reference_path,out_path):
    #Function to run make a nucl database and run blastn
    blastdb_name = reference_path + ".nhr"
    if not os.path.isfile(blastdb_name):
        print ("No blastdb found - running makeblastdb...")
        subprocess.call("makeblastdb -dbtype nucl -in {}".format(reference_path), shell=True)
    print("Running blastn alignment...")
    subprocess.call("blastn -outfmt 6 -perc_identity 95 -query {} -db {} > {}".format(query_path,reference_path,out_path), shell=True)
    print("Done!")
    return

def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def filter_top_blastn_hit(blastn_table):
    #function to produce a dictionary of dictionaries out of the top blastn hits
    print("Filtering blastn info...")
    blastn_dict = {}
    column_id_dict = {0:"query_id",1:"subject_id",2:"perc_id",3:"aln_len",\
        4:"mismatches",5:"gap_opens",6:"q_start",7:"q_end",8:"s_start",\
        9:"s_end",10:"evalue",11:"bit_score"}
    with open(blastn_table,"r") as in_table:
        for line in in_table:
            line_list = line.split()
            query_id = line_list[0]
            evalue = line_list[10]
            bit_score = line_list[11]
            #populate the empty dictionary
            if query_id not in blastn_dict:
                blastn_dict[query_id] = {}
                for index,info in enumerate(line_list[1:]):
                    blastn_dict[query_id][column_id_dict[index + 1]] = info
            #if the new evalue is lower, overwrite entry with new info
            #elif blastn_dict[query_id]['evalue'] > evalue:
            elif blastn_dict[query_id]['evalue'] > evalue:
                blastn_dict[query_id] = {}
                for index,info in enumerate(line_list[1:]):
                    blastn_dict[query_id][column_id_dict[index + 1]] = info
    print ("Done!")
    return blastn_dict

#Run blastn
run_blastn(de_novo_assembly,reference_assembly,blastn_output)

#Maybe this non-redundant procedure misses repeat regions..?
blastn_dict = filter_top_blastn_hit(blastn_output)

#1. Load pathway coordinates as a dict
pathway_table = pd.read_csv(args['pathway_table'], header=None, sep="\t")
pathway_dict = {}
total_aln_len = 0
print("Calculating fragmentation of pathways...")
print("Listing contigs that align to slm pathway...")
for count,pathway in enumerate(pathway_table[0]):
    pathway_min = int(pathway_table[1][count])
    pathway_max = int(pathway_table[2][count])
    #pathway_dict[pathway] = pathway_min + "-" + pathway_max
    pathway_dict[pathway] = {}
    pathway_dict[pathway]['min'] = pathway_min
    pathway_dict[pathway]['max'] = pathway_max
    pathway_dict[pathway]['len'] = pathway_max - pathway_min #+ 1
    pathway_dict[pathway]['contig_list'] = []
    pathway_dict[pathway]['len_recovered'] = 0
    #blastn_dict made non-redundant by lowest evalue
    for query in blastn_dict:
        aln_len = int(blastn_dict[query]['aln_len'])
        #for the first iteration, tally up total aln length
        if count == 0: total_aln_len += aln_len
        s_start,s_end = int(blastn_dict[query]['s_start']),int(blastn_dict[query]['s_end'])
        #I think the logic must be breaking down here..
        #What if the best hit of a long contig that contains a pathway(s), but best alignment is outside of pathway coordinate range...
        #Might need to inverse this logic...
        #If the subject start or stop is in the range of the pathway coordinates
        #if s_start in range(pathway_min,pathway_max) or s_end in range(pathway_min,pathway_max):
        s_min,s_max = min([s_start,s_end]),max([s_start,s_end])
        if pathway_min in range(s_min,s_max) or pathway_max in range(s_min,s_max):
            #Add the query contig to the list
            pathway_dict[pathway]['contig_list'].append(query)
            #Show me the contig list for the slm pathway
            if pathway == "slm": print query
            #The overlap function, expects poth sets of min and max values
            #s_min,s_max = min([s_start,s_end]),max([s_start,s_end])
            overlap = getOverlap([pathway_min,pathway_max],[s_min,s_max])
            pathway_dict[pathway]['len_recovered'] += overlap
        #Use the old logic as well
        elif s_start in range(pathway_min,pathway_max) or s_end in range(pathway_min,pathway_max):
            pathway_dict[pathway]['contig_list'].append(query)
            #Show me the contig list for the slm pathway
            if pathway == "slm": print query
            #The overlap function, expects poth sets of min and max values
            #s_min,s_max = min([s_start,s_end]),max([s_start,s_end])
            overlap = getOverlap([pathway_min,pathway_max],[s_min,s_max])
            pathway_dict[pathway]['len_recovered'] += overlap

#3. Report alignment info
output_name = os.path.splitext(os.path.basename(de_novo_assembly))[0] + ".out"
with open(output_name, "w") as outfile:
    outfile.write("pathway\tnum_contigs\tlen_recovered\ttotal_len\tbreak_len\n")
    for pathway in pathway_dict:
        num_contigs = str(len(pathway_dict[pathway]['contig_list']))
        len_recovered = str(pathway_dict[pathway]['len_recovered'])
        total_len = str(pathway_dict[pathway]['len'])
        break_len = str(int(total_len) - int(len_recovered))
        outfile.write("\t".join([pathway,num_contigs,len_recovered,total_len,break_len]) + "\n")

print("Your output is in: {}".format(output_name))
#print("Total aln length of this de novo assembly to the reference genome is: {}".format(total_aln_len))
