===========
CompleConta
===========

Estimation of completeness/contamination with simple taxonomic classification of bacterial genomes (bins)

Introduction:
-------------

**Comple** teness and **Conta** mination estimation using EggNOG-profiles

This tool was developed for `PhenDB <http://phendb.org/>`_, where it provides a quick and rough estimation of both completeness, contamination and heterogeneity (similar to *CheckM*) and inference of taxonomic information that can be gathered from those marker genes (like in *Amphora*). The analysis includes a rudimental lowest common ancestor (LCA) approach with 90% majority rule. Currently the databases are build from bactNOG data of the EggNOG 4.5 release and taxonomy level 2 of EggNOG 5 (Bacteria). Please note, that the taxonomic tree only contains **bacteria**, which were used to create the EggNOG database. To run CompleConta you need to perform a classification with the bactNOG database, which might be time consuming. We suggest to use HMMER 3.1b2 or 3.2.1 for that purpose.

Dependencies:
-------------

* `BioPython <https://biopython.org/wiki/Download>`_ tested with v1.71
* `MUSCLE Aligner <https://www.drive5.com/muscle/>`_ tested with v3.8.31
* `ncbi-blast+ <https://blast.ncbi.nlm.nih.gov/Blast.cgi>`_ v2.3.0 or higher

Usage:
------
.. code-block::

    usage: compleconta.py [-h] [--margin MARGIN] [--majority MAJORITY]
                          [--rank RANK] [--aai AAI] [-o TAXONOMY_OUTPUT]
                          [--threads N_BLAST_THREADS] [--muscle MUSCLE_EXECUTABLE]
                          [--blast BLAST_EXECUTABLE] [--database DATABASE]
                          input.faa input.faa.out

    Completeness and Contamination estimation using EggNOG-profiles

    positional arguments:
      input.faa             genome proteome file
      input.faa.out         tab separated file with bactNOG classification of
                            proteome input file

    optional arguments:
      -h, --help            show this help message and exit
      --margin MARGIN       Taxonomy: fraction margin for hits taken relative to
                            bitscore of best hit in blastp (default: 0.9)
      --majority MAJORITY   Taxonomy: majority threshold for fraction of paths
                            that support the lowest common ancestor (default: 0.9)
      --rank RANK           Taxonomy: lowest standard rank to be reported. 0:
                            species, 1: genus, ... 5: phylum (default: 1)
      --aai AAI             Contamination: amino acid identity to which the
                            multiple marker genes are considered to be strain
                            heterogenic (default: 0.9)
      -o TAXONOMY_OUTPUT    Taxonomy: file to write additional taxonomic
                            information for each marker gene, if not set,
                            information is omitted (default: None)
      --threads N_BLAST_THREADS
                            Taxonomy: number of parallel blastp jobs run (default:
                            5)
      --muscle MUSCLE_EXECUTABLE
                            Path to the muscle executable (default: None)
      --blast BLAST_EXECUTABLE
                            Path to the blast executable (makeblastdb) (default:
                            None)
      --database DATABASE   database which was used for annotation (default: auto)


Positional argument 1 should be the fasta format proteome file, positional argument 2 the corresponding classification output as **tab-separated** file where column 1 is the sequence identifier of the fasta file and column 2 the ENOG cluster id to which it was classified. Example files can be found in the example directory.



Output:
-------

The output is written to :code:`stdout` and provides Completeness, Contamination, Strain Heterogeneity (all in fractions), ncbi taxon ID, taxon name and taxon rank.

The output of the example for EggNOG 5.0 is\*:

.. code-block::

    Comp.	Cont.	St. Het.	ncbi_taxid	taxon_name	taxon_rank
    1.0000	0.0000	0.0000	1760	Actinobacteria	class

\*Note that the result may differ with other classification tools or databases used.

The optional output file contains detailed information about the marker genes that were used and their taxonomic information. Note that the majority cutoff for proceding to a lower level is 0.9, the ranks below are still reported in the output file, but not used for the final LCA.

.. code-block::

    LCA path and percentage of marker genes assignment:
    Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 0.88	Nocardioidaceae 0.88

    LCA per sequence of identified marker genes:i
    WP_008360965.1	COG0048	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 0.99	Streptomycetales 0.31	Streptomycetaceae 0.31	Streptomyces 0.29
    WP_008360994.1	COG0051	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Streptomycetales 0.30	Streptomycetaceae 0.30	Streptomyces 0.27
    WP_008362894.1	COG0080	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 0.93	Nocardioidaceae 0.93	Nocardioides 0.36
    WP_008362892.1	COG0081	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00	Pimelobacter 0.50
    WP_008360995.1	COG0087	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00	Nocardioides 0.50
    WP_008360999.1	COG0089	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00	Nocardioides 0.50
    WP_008361002.1	COG0090	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00
    WP_008361005.1	COG0091	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00
    WP_008361007.1	COG0092	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00
    WP_008361015.1	COG0093	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Streptomycetales 0.33	Streptomycetaceae 0.33	Streptomyces 0.30
    WP_008361018.1	COG0094	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00	Nocardioides 0.50
    WP_040755997.1	COG0096	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00	Nocardioides 0.67
    WP_008361023.1	COG0097	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00	Pimelobacter 0.33
    WP_008361027.1	COG0098	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00
    WP_008363218.1	COG0130	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00
    WP_008363208.1	COG0184	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00	Nocardioides 0.50
    WP_008361004.1	COG0185	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Micrococcales 0.44	Intrasporangiaceae 0.35	Nocardioides 0.15
    WP_008361013.1	COG0186	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00
    WP_008361017.1	COG0198	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00
    WP_008361030.1	COG0200	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00
    WP_008361063.1	COG0203	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00
    WP_008363328.1	COG0228	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00
    WP_008362879.1	COG0244	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00
    WP_008361025.1	COG0256	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00
    WP_040755263.1	COG0261	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00
    WP_008363690.1	COG0268	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00
    WP_008359372.1	COG0292	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00	Nocardioides 0.50
    WP_008363321.1	COG0335	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00	Nocardioides 0.50
    WP_008355442.1	COG0359	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00
    WP_008355448.1	COG0360	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00
    WP_008358493.1	COG0536	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00
    WP_008358688.1	COG0691	Bacteria 1.00	Actinobacteria 1.00	Actinobacteria 1.00	Propionibacteriales 1.00	Nocardioidaceae 1.00

Setup:
------

At the moment, there is no special setup required.

.. code-block:: bash

    # Clone repository
    git clone https://github.com/phyden/compleconta
    # Change to directory
    cd compleconta
    # Run to display useage
    ./compleconta.py -h

Both the reduced taxonomy files (:code:`names.dmp` and :code:`nodes.dmp`) and the databases which were created from the bactNOG raw alignments are located in the data folder. The tool is ready to run, and will create the indices for the database files on execution if non existent. The script to prepare the database from EggNOG 4.5 is provided: :code:`prepare_blast_database.sh`. To include other databases this script requires slight adaptions.

Please file an issue or contact the author if you need assistance.
