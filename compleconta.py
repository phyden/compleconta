#!/usr/bin/env python

import os, sys
from compleconta import FileIO, Annotation, EnogLists, aminoAcidIdentity, Check, MarkerGeneBlast, ncbiTaxonomyTree

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

data_dir=IOobj.get_data_dir()

database_dir=data_dir+"/databases"

taxid_list=MarkerGeneBlast.getTaxidsFromSequences(database_dir,gc_subset)

taxonomy_dir=data_dir+"/taxonomy"

tree=ncbiTaxonomyTree.NcbiTaxonomyTree(taxonomy_dir)

lca_taxid, lca_name, lca_rank=tree.getLCA(taxid_list,rank=0,majority_threshold=0.9) #standard rank used: 0 (species), majority threshold 0.9


#result is a tuple containing (completeness(fraction), contamination(fraction))
print("Comp.\tCont.\tSt. Het.\tncbi_taxid\ttaxon_name\ttaxon_rank\n%.4f\t%.4f\t%.4f\t%i\t%s\t%s" %(float(result[0]),float(result[1]), aai, lca_taxid, lca_name, lca_rank))

