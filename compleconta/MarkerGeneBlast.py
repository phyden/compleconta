#!/usr/env python

import sys, os, tempfile, subprocess
import multiprocessing


def prepareFiles(gc, enog_list):
    """ function which creates a temporary directory and input, output filepairs for each sequence with the info to which enog it belongs """

    tmp_dir = tempfile.mkdtemp()

    outfiles = []
    inputfiles = []
    enogs_used = []
    seqs_used = []

    for enog in enog_list:
        seqs = gc.get_sequences_by_enog(enog)
        for seqid in seqs.keys():
            tmpfile_handler, tempfile_output = tempfile.mkstemp(dir=tmp_dir)
            outfiles.append(tempfile_output)

            tmpfile_handler, tempfile_input = tempfile.mkstemp(dir=tmp_dir)
            inputfiles.append(tempfile_input)

            enogs_used.append(enog)
            seqs_used.append(seqid)

            with open(tempfile_input, "w") as tmpfile_handler:
                tmpfile_handler.write(">%s\n%s\n" % (seqid, seqs[seqid]))

    return tmp_dir, outfiles, inputfiles, enogs_used, seqs_used


def run_blast_job(parameter_set):
    """ run function which is called by the multiprocessing pool """

    database, inputfile, outputfile, margin, blast_executable, makeblastdb_executable = parameter_set

    if check_database(database, makeblastdb_executable) == 0:
        exit_status = subprocess.call([blast_executable, "-db", database, "-query", inputfile, "-out", outputfile,
                                       "-outfmt", "6"])
        best_hit = readOutput(outputfile, margin)
        return best_hit

    else:
        return []


def check_database(database, blast_executable):
    """ checks existance of database file & creates index files in necessary """

    endings = ["phr", "pin", "psq"]

    try:
        time_db = os.path.getctime(database)
    except OSError:
        sys.stderr.write("Error: database file not existing or inaccessable: %s\n" % database)
        return 1

    recreate = False

    for ending in endings:
        indexfile = ".".join([database, ending])
        try:
            if time_db > os.path.getctime(indexfile):
                recreate = True
                break
        except OSError:
            if os.path.isfile(indexfile):
                sys.stderr.write("Error: database index existing but inaccessable: %s\n" % indexfile)
                return 1
            else:
                recreate = True

    if recreate == True:
        sys.stderr.write("Info: database indices will be created for %s\n" % database)
        subprocess.call([blast_executable, "-in", database, "-dbtype", "prot"])

    return 0


def readOutput(outputfile, margin):
    """ function to parse the blastp tabular output and return the best hits (with a margin from the best bitscore downwards) """

    with open(outputfile, "r") as tmpfile_handler:
        maxscore = 0
        tophit = []
        for line in tmpfile_handler:
            # taxid, str_score = line.strip().split("\t")
            fields = line.strip().split("\t")
            hit = fields[1]
            pident = float(fields[2])
            bitscore = float(fields[11])
            maxscore = max(maxscore, bitscore)
            if bitscore < maxscore * margin:
                break
            else:
                tophit.append(hit)
        if len(tophit) == 0:
            tophit.append(str(1))
    return tophit


def getTaxidsFromSequences(databasepath, gc, args, blast_executable, makeblastdb_executable):
    """ master function that runs above functions parallelized in a worker pool of n_blast_threads subprocesses """

    n_blast_threads = max(args.n_blast_threads, 1)  # minimum number of workers: 1
    margin = min(abs(args.margin), 1.0)  # margin has to be in range 0.0-1.0

    best_hits = []

    enog_list = gc.get_profile()

    pool = multiprocessing.Pool(n_blast_threads)

    tmp_dir, outfiles, inputfiles, enog_list, seq_list = prepareFiles(gc, enog_list)

    parameter_sets = []

    for i in range(0, len(enog_list)):
        database = databasepath + "/" + enog_list[i] + ".fa"
        parameter_sets.append((database, inputfiles[i], outfiles[i], margin, blast_executable, makeblastdb_executable))
        # best_hit=runBlastJob(parameter_sets[i])

    best_hits = pool.map(run_blast_job, parameter_sets)

    os.system("rm -r %s" % tmp_dir)

    return best_hits, seq_list, enog_list
