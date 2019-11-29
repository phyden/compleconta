#!/usr/bin/env python

import os
import sys
from Bio import SeqIO

def load_sequences(protein_file):

    seq_return={}

    if os.path.isfile(protein_file):
        with open(protein_file) as infile:
            for record in SeqIO.parse(infile,"fasta"):
                seq_return[record.id]=str(record.seq)

    if len(seq_return) == 0:
        sys.stderr.write("ERROR: provided protein.fasta file empty: {}\n".format(protein_file))
        raise EOFError
    return seq_return


def load_enog_annotation(hmmer_outfile):

    proteins={}

    if os.path.isfile(hmmer_outfile):
        with open(hmmer_outfile) as infile:
            for line in infile:
                line=line.strip().split("\t")
                proteins[line[0]]=line[1]

    if len(proteins) == 0:
        sys.stderr.write("ERROR: provided genotype file empty: {}\n".format(hmmer_outfile))
        raise EOFError
    return proteins


def check_database(arg_database, sample_enogs):
    """
    Function to detect wheter the database specified as parameter is existing.
    If 'auto' is specified, compleconta tries to automatically identify the database used
    :param arg_database: String as provided by user
    :param sample_enogs: list of enogs found in the sample
    :return: string of database folder
    """

    if arg_database == 'auto':
        database = determine_database(sample_enogs)

    else:
        if os.path.exists(os.path.dirname(__file__)+"/../data/" + arg_database):
            database = arg_database
        else:
            database = None

    return database


def determine_database(sample_enogs):
    """
    Function to automatically determine which of the databases might have been used
    :param sample_enogs: list of all enogs identifiers found in the sample
    :return: string of database folder
    """

    compare_set = set(sample_enogs)
    data_dir = os.path.dirname(__file__)+"/../data"
    databases = os.listdir(data_dir)
    max_matches = 0
    database = None
    for db in databases:
        with open(data_dir+"/"+db+"/set_of_enogs.txt") as inf_h:

            enog_set = set([enog.strip() for enog in inf_h.readlines()])
            matches = len(enog_set.intersection(compare_set))
            if matches > max_matches:
                database = db
                max_matches = matches

    if database is None:
        database = "eggnog5"
        sys.stderr.write("WARNING: database could not be determined, using eggnog5\n")
    else:
        sys.stderr.write("INFO: database automatically determined. Using {}\n".format(database))

    return database


class FileIO:

    def __init__(self, eggnog_version="eggnog4", output_filename="output_file.tsv"):

    #set the environment:
        self.proj_dir=os.path.dirname(__file__)

        self.data_dir=self.proj_dir+"/../data/"+eggnog_version

        self.universal_cogs_file=self.data_dir+"/set_of_enogs.txt"
        self.sorted_enogs_file=self.data_dir+"/copynumber_counts.tsv"

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

