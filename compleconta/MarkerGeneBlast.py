#!/usr/env python

import sys, os, tempfile
import multiprocessing

import Annotation

def prepareFiles(gc,enog_list):
    
    tmp_dir=tempfile.mkdtemp()
    
    outfiles=[]
    inputfiles=[]
    enogs_used=[]
    seqs_used=[]

    for enog in enog_list:
        seqs = gc.get_sequences_by_enog(enog)
        for seqid in seqs.keys():
            
            tmpfile_handler, tempfile_output = tempfile.mkstemp(dir=tmp_dir)
            outfiles.append(tempfile_output)

            tmpfile_handler, tempfile_input = tempfile.mkstemp(dir=tmp_dir)
            inputfiles.append(tempfile_input)

            enogs_used.append(enog)
	    seqs_used.append(seqid)

            with open(tempfile_input,"w") as tmpfile_handler:
                tmpfile_handler.write(">%s\n%s\n" % (seqid, seqs[seqid]))

    return tmp_dir, outfiles, inputfiles, enogs_used, seqs_used


def runBlastJob(parameter_set):

    database, inputfile, outputfile = parameter_set

    os.system("blastp -db %s -query %s -out %s -outfmt '6'" % parameter_set)

    best_hit = readOutput(outputfile)
    
    return best_hit

def readOutput(outputfile):

    with open(outputfile,"r") as tmpfile_handler:
        maxscore=0
        tophit=[]
        for line in tmpfile_handler:
            #taxid, str_score = line.strip().split("\t")
            fields=line.strip().split("\t")
            hit=fields[1]
            pident=float(fields[2])
	    bitscore=float(fields[11])
            maxscore=max(maxscore,bitscore)
            if bitscore < maxscore*0.9:
                break
            else:
                tophit.append(hit)
        if len(tophit)==0:
            tophit.append(str(1))
    return tophit


def getTaxidsFromSequences(databasepath,gc):

    best_hits=[]

    enog_list=gc.get_profile()

    pool=multiprocessing.Pool(5)

    tmp_dir, outfiles, inputfiles, enog_list, seq_list = prepareFiles(gc, enog_list)
    

    parameter_sets=[]
    for i in range(0,len(enog_list)):
        database=databasepath+"/"+enog_list[i]+".fa"
        parameter_sets.append((database, inputfiles[i], outfiles[i]))
        
    best_hits=pool.map(runBlastJob,parameter_sets)

    os.system("rm -r %s" % tmp_dir)

    return best_hits, seq_list, enog_list
 
