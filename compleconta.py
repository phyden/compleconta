#!/usr/bin/env python

import os, sys
from compleconta import FileIO, Annotation, EnogLists, aminoAcidIdentity, Check

# usage: compleconta.py /path/to/protein_file.faa /path/to/hmmer_results.faa.out

protein_file=sys.argv[1]
hmmer_file=sys.argv[2]

#using the FileIO class from compleconta. the enog lists+enog weights are stored in two files which are also found in compleconta/data
IOobj=FileIO.FileIO()

#function read_enog_list returns a list only if no header present (first column), or a dict additionally (all information)
all_enogs, enog_dict = IOobj.read_enog_list(IOobj.sorted_enogs_file,header=True)
curated34_list = IOobj.read_enog_list(IOobj.universal_cogs_file,header=False)

#to handle the weights, I created a EnogList class. initialized with the enog list and the dictionary
marker_set=EnogLists.EnogList(curated34_list, enog_dict)

#the genecollection contains all enogs, the sequence names associated and the sequences. assumption: inputfile = <proteins>.faa and hmmer classification results in <proteins>.faa.out, same directory
gc=Annotation.GeneCollection()
gc.create_from_file(protein_file, hmmer_file)

#subset to enogs that actually are in the list - needed for AAI, speeds up cc slightly
gc_subset=gc.subset(curated34_list)


aai=aminoAcidIdentity.aai_check(0.9,gc_subset)
result=Check.check_genome_cc_weighted(marker_set,gc.get_profile())

#result is a tuple containing (completeness(fraction), contamination(fraction))
print("Comp.\tCont.\tSt. Het.\n%.4f\t%.4f\t%.4f" %(float(result[0]),float(result[1]), aai))

