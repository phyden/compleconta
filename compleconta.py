#!/usr/bin/env python

import sys, subprocess
from compleconta import FileIO, Annotation, EnogLists, aminoAcidIdentity, Check, MarkerGeneBlast, ncbiTaxonomyTree

import argparse

parser = argparse.ArgumentParser(description='Completeness and Contamination estimation using EggNOG-profiles', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('protein_file', metavar='input.faa', type=str,
                    help='genome proteome file')
parser.add_argument('hmmer_file', metavar='input.faa.out', type=str,
                    help='tab separated file with bactNOG classification of proteome input file')
parser.add_argument('--margin', dest='margin', type=float, default=0.9,
                    help='Taxonomy: fraction margin for hits taken relative to bitscore of best hit in blastp')
parser.add_argument('--majority', dest='majority', type=float, default=0.9,
                    help='Taxonomy: majority threshold for fraction of paths that support the lowest common ancestor')
parser.add_argument('--rank', dest='rank', type=int, default=1,
                    help='Taxonomy: lowest standard rank to be reported. 0: species, 1: genus, ... 5: phylum')
parser.add_argument('--aai', dest='aai', type=float, default=0.9,
                    help='Contamination: amino acid identity to which the multiple marker genes are considered to be strain heterogenic')
parser.add_argument('-o', dest='taxonomy_output', type=str, required=False,
                    help='Taxonomy: file to write additional taxonomic information for each marker gene, if not set, information is omitted')
parser.add_argument('--threads', dest='n_blastp_threads', type=int, default=5,
                    help='Taxonomy: number of parallel blastp jobs run')


args = parser.parse_args()

# usage: compleconta.py /path/to/protein_file.faa /path/to/hmmer_results.faa.out

required_executables=["blastp","makeblastdb","muscle"]
for requirement in required_executables:
    try:
        status=subprocess.check_output([requirement,"-version"])
    except OSError:
        sys.stderr.write("Error: executable '%s' not found in PATH\nAborting\n" % requirement)
        exit(1)
        

protein_file=args.protein_file
hmmer_file=args.hmmer_file

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


aai=aminoAcidIdentity.aai_check(args.aai,gc_subset)
completeness, contamination=Check.check_genome_cc_weighted(marker_set,gc.get_profile())

data_dir=IOobj.get_data_dir()

database_dir=data_dir+"/databases"

taxid_list, sequence_ids, enog_names=MarkerGeneBlast.getTaxidsFromSequences(database_dir,gc_subset, args.n_blastp_threads, args.margin)

taxonomy_dir=data_dir+"/taxonomy"

tree=ncbiTaxonomyTree.NcbiTaxonomyTree(taxonomy_dir)

lca_per_sequence=[]
nodes_per_sequence=[]
percentages_per_sequence=[]
for sub_taxids in taxid_list:
	reported_lca, nodes, percentages=tree.getLCA(sub_taxids,rank=args.rank,majority_threshold=args.majority)
	lca_per_sequence.append(reported_lca.taxid)
	nodes_per_sequence.append(nodes)
	percentages_per_sequence.append(percentages)
		

reported_lca, nodes, percentages=tree.getLCA(lca_per_sequence,rank=args.rank,majority_threshold=args.majority) #standard ranks: 0 (species), 1 (genus), ..., majority threshold 0.9


#result is a tuple containing (completeness(fraction), contamination(fraction))
print("Comp.\tCont.\tSt. Het.\tncbi_taxid\ttaxon_name\ttaxon_rank\n%.4f\t%.4f\t%.4f\t%i\t%s\t%s" %(float(completeness),float(contamination), aai, reported_lca.taxid, reported_lca.name, reported_lca.rank))

if args.taxonomy_output:
    output_file=args.taxonomy_output
    with open(output_file, "w") as outfile_handler:
        outfile_handler.write("LCA path and percentage of marker genes assignment:\n")
        outfile_handler.write("%s\n" % "\t".join([nodes[i].name+" "+str(round(percentages[i],2)) for i in range(len(nodes))]))

        outfile_handler.write("\nLCA per sequence of identified marker genes:i\n")
        for i in range(len(sequence_ids)):
	    taxonomy="\t".join([nodes_per_sequence[i][j].name+" "+str(round(percentages_per_sequence[i][j],2)) for j in range(len(nodes_per_sequence[i]))])
	    outfile_handler.write("%s\t%s\t%s\n"%(sequence_ids[i],enog_names[i],taxonomy))
