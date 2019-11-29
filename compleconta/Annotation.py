#!/usr/env python

from compleconta import FileIO


class GeneCollection:

    def __init__(self):

        self.enogs = []
        self.enog_to_genes = {}
        self.genes_to_enog = {}
        self.sequences = {}

    def create_from_file(self, protein_file, hmmer_outfile, genome_id="NA"):
        self.id = genome_id

        self.load_sequences(protein_file)
        self.load_enog_annotation(hmmer_outfile)

        self.enogs = []
        self.enog_to_genes = {}

        for gene in self.genes_to_enog.keys():
            self.enogs.append(self.genes_to_enog[gene])
            if not self.enog_to_genes.get(self.genes_to_enog[gene]):
                self.enog_to_genes[self.genes_to_enog[gene]] = []
            self.enog_to_genes[self.genes_to_enog[gene]].append(gene)

    def subset(self, enogs=None):

        # return a subset of GeneCollection object and return a new one without changing the original set

        new_set = GeneCollection()
        new_set.id = self.id

        if not enogs:
            enogs = self.enogs

        new_set.enogs = enogs[:]

        for enog in enogs:
            if self.enog_to_genes.get(enog):
                new_set.enog_to_genes[enog] = self.enog_to_genes.get(enog)
                for gene in new_set.enog_to_genes.get(enog, []):
                    new_set.genes_to_enog[gene] = enog
                    new_set.sequences[gene] = self.sequences.get(gene)

        return new_set

    def get_profile(self):

        return self.enogs

    def get_sequences_by_enog(self, enog):

        seqs = {}
        for gene in self.enog_to_genes.get(enog, []):
            seqs[gene] = self.sequences[gene]

        return seqs

    def load_sequences(self, protein_file):

        self.sequences = FileIO.load_sequences(protein_file)

    def load_enog_annotation(self, readfile):

        self.genes_to_enog = FileIO.load_enog_annotation(readfile)

    def get_multicopy_enogs(self):

        mc_enogs = []
        for enog in self.enogs:
            if len(self.get_sequences_by_enog(enog)) > 1:
                mc_enogs.append(enog)

        return mc_enogs


class Gene:

    def __init__(self, gene_id, enog, sequence):
        self.id = gene_id
        self.enog = enog
        self.sequence = sequence

    def get_sequence():
        return self.sequence

    def get_enog():
        return self.enog
