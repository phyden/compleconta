###############################################################################
#
# aminoAcidIdentity.py - calculate AAI between aligned marker genes
# copied from checkM. modified by Patrick Hyden, Oct 2017
#
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import os
import sys
import subprocess
from Bio import AlignIO

def make_alignments(sequences):

    tmpfasta=[]
    for header in sequences.keys():
        tmpfasta.append(">"+header)
        tmpfasta.append(sequences[header])

    muscle_executable="/apps/muscle/3.8.31/muscle"
    child=subprocess.Popen(muscle_executable,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True,shell=False)
    child.stdin.write("\n".join(tmpfasta))
    child.stdin.close()
    alignment = AlignIO.read(child.stdout, "fasta")

    return str(alignment[0].seq), str(alignment[1].seq)


def aai_check(aaiStrainThreshold, gene_collection):
    """Calculate AAI between input alignments."""
    aaiRawScores = {}
    aaiHetero = {}
    aaiMeanBinHetero = {}

    mc_enogs=gene_collection.get_multicopy_enogs()
    for markerId in mc_enogs:
        seqs=gene_collection.get_sequences_by_enog(markerId)        
        aaiRawScores[markerId]=[]
        for i in range(len(seqs)):
            seqIdI, seqI = seqs.popitem()
            for seqIdJ, seqJ in seqs.items():
                seqI, seqJ = make_alignments({seqIdI: seqI, seqIdJ: seqJ})
                aai = aai_seq(seqI, seqJ)
                aaiRawScores[markerId].append(aai)

    aaiHetero, aaiMeanBinHetero = strainHetero(aaiRawScores, aaiStrainThreshold)

    return aaiMeanBinHetero

def strainHetero(aaiScores, aaiStrainThreshold):
    """Calculate strain heterogeneity."""
    aaiMeanBinHetero = {}
    aaiHetero = {}
    strainCount = 0
    multiCopyPairs = 0


    for markerId in aaiScores.keys():
        localStrainCount = 0
        for aaiScore in aaiScores[markerId]:
            multiCopyPairs += 1
            if aaiScore > aaiStrainThreshold:
                strainCount += 1
                localStrainCount += 1
            strainHetero = float(localStrainCount) / len(aaiScores[markerId])
            aaiHetero[markerId] = strainHetero

    if not multiCopyPairs==0:
        aaiMeanBinHetero = float(strainCount) / multiCopyPairs
    else:
        aaiMeanBinHetero = 0

    return aaiHetero, aaiMeanBinHetero

def aai_seq(seq1, seq2):
    """Calculate amino acid identity between sequences."""
    assert len(seq1) == len(seq2)

    # calculation of AAI should ignore missing data at
    # the start of end of each sequence
    startIndex = 0
    for i in range(0, len(seq1)):
        if seq1[i] == '-' or seq2[i] == '-':
            startIndex = i + 1
        else:
            break

    endIndex = len(seq1)
    for i in range(len(seq1) - 1, 0, -1):
        if seq1[i] == '-' or seq2[i] == '-':
            endIndex = i
        else:
            break

    mismatches = 0
    seqLen = 0
    for i in range(startIndex, endIndex):
        if seq1[i] != seq2[i]:
            mismatches += 1
            seqLen += 1
        elif seq1[i] == '-' and seq2[i] == '-':
            pass
        else:
            seqLen += 1

    if seqLen == 0:
        return 0.0

    return 1.0 - (float(mismatches) / seqLen)
