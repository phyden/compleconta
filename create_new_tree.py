#!/usr/bin/env python2

import tempfile
import subprocess

from compleconta.ncbiTaxonomyTree import NcbiTaxonomyTree
from compleconta.FileIO import FileIO

ncbi_taxonomy_tar = "/mirror/ncbi/current/taxonomy/taxdump.tar.gz"
tmp_dir = tempfile.mkdtemp()

subprocess.call(["tar", "-xf", ncbi_taxonomy_tar, "-C", tmp_dir])

environment = FileIO()
data_dir = environment.get_data_dir()
inputPathTaxonomy = tmp_dir

outputPathTaxonomy = data_dir + "/taxonomy"

selectionFile = data_dir + "/tax_ids_used.txt"

tree = NcbiTaxonomyTree(inputPathTaxonomy)

selected_nodes = []

with open(selectionFile, "r") as inf_h:
    for line in inf_h:
        taxid = int(line.strip())
        try:
            subtree = tree.getAscendantsWithRanksAndNames([taxid], only_std_ranks=False)  # [taxid]

            for node in subtree[taxid]:
                selected_nodes.append(node.taxid)
        except KeyError:
            print(taxid)

lookup = set(selected_nodes)

newNodesFileName = outputPathTaxonomy + "/nodes.dmp"
newNodesFile = open(newNodesFileName, "w")

with open(inputPathTaxonomy + "/nodes.dmp") as inf_h:
    for line in inf_h:
        fields = line.split("|")
        taxid = int(fields[0].strip())
        if taxid in lookup:
            newNodesFile.write(line)

newNodesFile.close()

newNamesFileName = outputPathTaxonomy + "/names.dmp"
newNamesFile = open(newNamesFileName, "w")

with open(inputPathTaxonomy + "/names.dmp") as inf_h:
    for line in inf_h:
        fields = line.split("|")
        taxid = int(fields[0].strip())
        name_type = fields[3].strip()
        if taxid in lookup and name_type == "scientific name":
            newNamesFile.write(line)

newNamesFile.close()

subprocess.call(["rm", "-r", tmp_dir])
