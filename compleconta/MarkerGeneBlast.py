#!/usr/env python

import sys, os, tempfile
import Annotation

def getTaxidsFromSequences(databasepath,gc):

    #load module ncbiblastplus (blastp must be in path!)
    os.system("module load ncbiblastplus")

    best_hits=[]

    tmp_dir=tempfile.mkdtemp()

    outfiles=[]
    
    enog_list=gc.get_profile()
    for enog in enog_list:
	#with open(tempfile_input,"w") as tmpfile_handler:
	tmpfile_handler, tempfile_output = tempfile.mkstemp(dir=tmp_dir)
        outfiles.append(tempfile_output)

        tmpfile_handler, tempfile_input = tempfile.mkstemp(dir=tmp_dir)

        with open(tempfile_input,"w") as tmpfile_handler:
            seqs = gc.get_sequences_by_enog(enog)
            for seqid in seqs.keys():
                tmpfile_handler.write(">%s\n%s\n" % (seqid, seqs[seqid]))

	database=databasepath+"/"+enog+".fa"
        os.system("blastp -db %s -query %s -out %s -outfmt '6 sseqid score'" % (database, tempfile_input, tempfile_output))
    

    for tempfile_output in outfiles:
	
	with open(tempfile_output,"r") as tmpfile_handler:
            maxscore=0
            tophit=""
            for line in tmpfile_handler:
                taxid, str_score = line.strip().split("\t")
                score=float(str_score)
                if score > maxscore:
                    maxscore=score
                    tophit=taxid

        best_hits.append(tophit)
    os.system("rm -r %s" % tmp_dir)
    
    return best_hits 
