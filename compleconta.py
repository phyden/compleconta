#!/usr/bin/env python3

# usage: compleconta.py /path/to/protein_file.faa /path/to/hmmer_results.faa.out

import sys
import os
import subprocess
import argparse

from compleconta import FileIO, Annotation, EnogLists, aminoAcidIdentity, Check, MarkerGeneBlast, ncbiTaxonomyTree


def get_args():
    parser = argparse.ArgumentParser(description='Completeness and Contamination estimation using EggNOG-profiles',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

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
    parser.add_argument('--threads', dest='n_blast_threads', type=int, default=5,
                        help='Taxonomy: number of parallel blastp jobs run')
    parser.add_argument('--muscle', dest='muscle_executable', type=str, required=False,
                        help='Path to the muscle executable')
    parser.add_argument('--blast', dest='blast_executable', type=str, required=False,
                        help='Path to the blast executable (makeblastdb)')

    # check for required executables

    return parser.parse_args()


def check_requirements(args):
    """Simple function that checks for BLAST and MUSCLE executables"""

    if args.blast_executable:
        if not os.path.exists(args.blast_executable):
            blastp = "blast"
        else:
            if os.path.isdir(args.blast_executable):
                blastp = os.path.join(args.blast_executable, "blastp")
                makeblastdb = os.path.join(args.blast_executable, "makeblastdb")
            else:
                if os.path.basename(args.blast_executable) == "blastp":
                    blastp = args.blast_executable
                    makeblastdb = os.path.join(os.path.dirname(args.blast_executable),"makeblastdb")
                elif os.path.basename(args.blast_executable) == "makeblastdb":
                    blastp = os.path.join(os.path.dirname(args.blast_executable), "blastp")
                    makeblastdb = args.blast_executable

    if args.muscle_executable:
        if not os.path.exists(args.muscle_executable):
            muscle = "muscle"
        else:
            if os.path.isdir(args.muscle_executable):
                muscle = os.path.join(args.muscle_executable, muscle)
            else:
                muscle = args.muscle_executable

    required_executables = (blastp, makeblastdb, muscle)

    for requirement in required_executables:
        try:
            status = subprocess.check_output([requirement, "-version"])
        except OSError:
            sys.stderr.write('Error: executable \'{module}\' not found in PATH\nExiting\n'.format(module=requirement))
            exit(1)

    return required_executables


def main():
    """Main function"""

    args = get_args()
    blast_executable, makeblastdb_executable, muscle_executable = check_requirements(args)

    protein_file = args.protein_file
    hmmer_file = args.hmmer_file

    # using the FileIO class from compleconta. the enog lists+enog weights are stored in two files which are also found in compleconta/data
    IOobj = FileIO.FileIO()

    # function read_enog_list returns a list only if no header present (first column), or a dict additionally (all information)
    all_enogs, enog_dict = IOobj.read_enog_list(IOobj.sorted_enogs_file, header=True)
    curated34_list = IOobj.read_enog_list(IOobj.universal_cogs_file, header=False)

    # to handle the weights, I created a EnogList class. initialized with the enog list and the dictionary
    marker_set = EnogLists.EnogList(curated34_list, enog_dict)

    # the genecollection contains all enogs, the sequence names associated and the sequences. assumption: inputfile = <proteins>.faa and hmmer classification results in <proteins>.faa.out, same directory
    gc = Annotation.GeneCollection()
    gc.create_from_file(protein_file, hmmer_file)

    # subset to enogs that actually are in the list - needed for AAI, speeds up cc slightly
    gc_subset = gc.subset(curated34_list)

    aai = aminoAcidIdentity.aai_check(gc_subset, args, muscle_executable)
    completeness, contamination = Check.check_genome_cc_weighted(marker_set, gc.get_profile())

    data_dir = IOobj.get_data_dir()

    database_dir = data_dir + "/databases"

    taxid_list, sequence_ids, enog_names = MarkerGeneBlast.get_taxids_of_sequences(database_dir, gc_subset, args,
                                                                                   blast_executable, makeblastdb_executable)

    taxonomy_dir = data_dir + "/taxonomy"

    tree = ncbiTaxonomyTree.NcbiTaxonomyTree(taxonomy_dir)

    lca_per_sequence = []
    nodes_per_sequence = []
    percentages_per_sequence = []
    for sub_taxids in taxid_list:
        reported_lca, nodes, percentages = tree.getLCA(sub_taxids, rank=args.rank, majority_threshold=args.majority)
        lca_per_sequence.append(reported_lca.taxid)
        nodes_per_sequence.append(nodes)
        percentages_per_sequence.append(percentages)

    # standard ranks: 0 (species), 1 (genus), ..., majority threshold 0.9
    reported_lca, nodes, percentages = tree.getLCA(lca_per_sequence, rank=args.rank, majority_threshold=args.majority)

    # result is a tuple containing (completeness(fraction), contamination(fraction))
    sys.stdout.write("Comp.\tCont.\tSt. Het.\tncbi_taxid\ttaxon_name\ttaxon_rank\n")
    sys.stdout.write("{:.4f}\t{:.4f}\t{:.4f}\t{}\t{}\t{}\n".format(float(completeness), float(contamination), aai,
                                                                reported_lca.taxid, reported_lca.name,
                                                                reported_lca.rank))

    if args.taxonomy_output:
        output_file = args.taxonomy_output
        with open(output_file, "w") as outfile_handler:
            outfile_handler.write("LCA path and percentage of marker genes assignment:\n")
            outfile_handler.write("\t".join(["{} {:.2f}".format(n.name, p) for n, p in zip(nodes, percentages)]))

            outfile_handler.write("\n\nLCA per sequence of identified marker genes:i\n")
            for i in range(len(sequence_ids)):
                taxonomy = "\t".join(["{} {:.2f}".format(node.name, perc) for node, perc in
                                      zip(nodes_per_sequence[i], percentages_per_sequence[i])])
                outfile_handler.write("{}\t{}\t{}\n".format(sequence_ids[i], enog_names[i], taxonomy))


if __name__ == "__main__":
    main()
