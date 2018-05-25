# CompleConta
**Comple**teness and **Conta**mination estimation using EggNOG-profiles

This tool was developed for [PhenDB](http://phendb.org/), where it provides a quick and rough estimation of both completeness, contamination and heterogeneity (similar to *CheckM*) and inference of taxonomic information that can be gathered from those marker genes (like in *Amphora*). The analysis includes a rudimental lowest common ancestor (LCA) approach with 90% majority rule. Currently the databases are build from bactNOG data of the EggNOG 4.5 release. Please note, that the taxonomic tree only contains nodes present in bactNOG. To run CompleConta you need to perform a classification with the bactNOG database, which might be time consuming. We suggest to use HMMER 3.1 for that purpose.

## Dependencies:

* [BioPython](https://biopython.org/wiki/Download) tested with v1.71
* [MUSCLE Aligner](https://www.drive5.com/muscle/) tested with v3.8.31
* [ncbi-blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi) v2.3.0 or higher

## Usage:

`compleconta.py input.faa classification.csv`

Positional argument 1 should be the fasta format proteome file, positional argument 2 the corresponding classification output as <tab-separated> file where column 1 is the sequence identifier of the fasta file and column 2 the ENOG cluster id to which it was classified. Example files can be found in the example directory.

## Output:

The output is written to `stdout` and provides Completeness, Contamination, Strain Heterogeneity (all in fractions), ncbi taxon ID, taxon name and taxon rank. Tab separated with header:
```
Comp.	Cont.	St. Het.	ncbi_taxid	taxon_name	taxon_rank
1.0000	0.0294	0.0000	85015	Nocardioidaceae	family
```

The optional output file contains detailed information about the marker genes that were used and their taxonomic information. Note that the majority cutoff for proceding to a lower level is 0.9, the ranks below are still reported in the output file, but not used for the final LCA.

```
LCA path and percentage of marker genes assignment:
Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 0.91	Nocardioidaceae 0.91

LCA per sequence of identified marker genes:
WP_008358950.1	ENOG4105C3G	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
WP_008359388.1	ENOG4105CSS	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
WP_008357500.1	ENOG4105CSS	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
WP_008360965.1	ENOG4108UKE	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 0.99	Corynebacteriales 0.27	Streptomycetaceae 0.25	Streptomyces 0.24
WP_008360967.1	ENOG4108UHY	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0	Nocardioides 0.5
WP_008363303.1	ENOG4105CE9	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
WP_008362894.1	ENOG4108UIK	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0	Kribbella 0.25
WP_008362892.1	ENOG4105C64	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
WP_008360995.1	ENOG4105EEE	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
WP_008360997.1	ENOG4106U5A	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
WP_008361002.1	ENOG4105CFD	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
WP_008361007.1	ENOG4105CKE	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
WP_008361015.1	ENOG4108UNN	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Micrococcales 0.32	Streptomycetaceae 0.25	Streptomyces 0.25
WP_008361018.1	ENOG4105CW6	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
WP_040755997.1	ENOG4108UJY	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0	Nocardioides 0.5
WP_008361023.1	ENOG4108R5J	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
WP_008361027.1	ENOG4108RA9	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
WP_008361056.1	ENOG4108Z04	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
WP_008361058.1	ENOG4108UHH	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
WP_008362447.1	ENOG4108UM5	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
WP_040756861.1	ENOG4108UJD	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
WP_008355939.1	ENOG4105CGR	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
WP_008363208.1	ENOG4105K77	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0	Nocardioides 0.5
WP_008361004.1	ENOG4105K7S	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 0.5	Nocardioidaceae 0.33	Janibacter 0.17
WP_008361013.1	ENOG4105K87	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
WP_008361009.1	ENOG4108R70	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0	Nocardioides 0.5
WP_008361030.1	ENOG4108UZ0	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
WP_008361031.1	ENOG4105CGG	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
WP_008361061.1	ENOG4105CTF	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0	Nocardioides 0.5
WP_008361025.1	ENOG4105K4C	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
WP_008363678.1	ENOG4105C8T	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
WP_008361060.1	ENOG4105G6W	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0	Nocardioides 0.5
WP_008358562.1	ENOG4105CA4	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
WP_008362431.1	ENOG4105CPM	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
WP_008363336.1	ENOG4105CB9	Bacteria 1.0	Actinobacteria 1.0	Actinobacteria 1.0	Propionibacteriales 1.0	Nocardioidaceae 1.0
```

## Setup:

Both the reduced taxonomy files (names.dmp and nodes.dmp) and the databases which were created from the bactNOG raw alignments are located in the data folder. The tool is ready to run, and will create the indices for the database files on execution if non existent. The indexing will also be done when `prepare_blast_database.sh` is run.

If you wish to expand or change to a different database, you need to specify a different set of marker genes. Since the ENOG-ids should be kept for newer versions of EggNOG-DB (whenever it will be released), it might be sufficient to specifiy the paths for the bactNOG alignment tar-ball and the taxonomy dump tar, and run the `prepare_blast_database.sh` script. If you run into troubles please open an issue or contact the author.
