#!/usr/bin/env python

import sys, os
from Bio import SeqIO
#import numpy as np
#import numpy.random as random

def load_sequences(protein_file):

    seq_return={}


    if os.path.isfile(protein_file):
        with open(protein_file) as infile:
            for record in SeqIO.parse(infile,"fasta"):
                seq_return[record.id]=str(record.seq)

    return seq_return


def load_enog_annotation(hmmer_outfile):

    proteins={}

    if os.path.isfile(hmmer_outfile):
        with open(hmmer_outfile) as infile:
            for line in infile:
                line=line.strip().split("\t")
                proteins[line[0]]=line[1]

    return proteins


class FileIO:

    def __init__(self, output_filename="output_file.tsv"):

    #set the environment:
        self.proj_dir=os.path.dirname(__file__)

        self.data_dir=self.proj_dir+"/../data"

        self.universal_cogs_file=self.data_dir+"/curated_34_enogs.txt"
        self.sorted_enogs_file=self.data_dir+"/sc_mc_counts_sorted_sc_short.tsv"

    ## functions ##

    def get_data_dir(self):
        return self.data_dir

    def read_enog_list(self, inputfname, header):
    
        # get sorted list of ENOGs ignoring header or without header
        read_list=[]
    
        if header==True:
            first_line=True
            read_dict={}
            header_list=[]
        else:
            first_line=False

        with open(inputfname,"r") as readfile:
            for line in readfile:
                if first_line==True:
                    first_line=False
                    header_list=line.strip().split("\t")
                else:
                    tmpline=line.strip().split("\t")
                    if header==True:
                        read_dict[tmpline[0]]={}
                        for i in range(1,len(tmpline)):
                           read_dict[tmpline[0]][header_list[i]]=tmpline[i]
                    #else:
                    read_list.append(tmpline[0])
    
        if header==True:
            return read_list, read_dict
        else:
            return read_list
    
    def write_results(self, results_dict):

        with open(self.output_file,"w") as outfile:
            outfile.write("%s\n" % "\t".join(results_dict.keys()))
            for i in range(len(results_dict[results_dict.keys()[0]])):
                outfile.write("%s\n" % "\t".join(str(results_dict[key][i]) for key in results_dict.keys()))

