###############################################################################
#
#   aminoAcidIdentity.py - calculate AAI between aligned marker genes
#   copied from checkM. modified by Patrick Hyden, Oct 2017
#   Made python3 compatible and substantial refactoring in Sept 2019
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

import subprocess
from Bio import AlignIO


def make_alignments(sequences, muscle_executable="muscle"):
    tmpfasta = []
    for header in sequences.keys():
        tmpfasta.append(">" + header)
        tmpfasta.append(sequences[header])

    child = subprocess.Popen(muscle_executable, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             universal_newlines=True, shell=False)
    child.stdin.write("\n".join(tmpfasta))
    child.stdin.close()
    alignment = AlignIO.read(child.stdout, "fasta")

    return str(alignment[0].seq), str(alignment[1].seq)


def aai_check(gene_collection, args):
    """Calculate AAI between input alignments."""
    aai_raw_scores = {}
    aai_strain_threshold = args.aai

    mc_enogs = gene_collection.get_multicopy_enogs()
    for markerId in mc_enogs:
        seqs = gene_collection.get_sequences_by_enog(markerId)
        aai_raw_scores[markerId] = []
        for i in range(len(seqs)):
            seq_id_i, seq_i = seqs.popitem()
            for seq_id_j, seq_j in seqs.items():
                seq_i, seq_j = make_alignments({seq_id_i: seq_i, seq_id_j: seq_j},
                                               muscle_executable=args.muscle_executable)
                aai = aai_seq(seq_i, seq_j)
                aai_raw_scores[markerId].append(aai)

    aai_hetero, aai_mean_bin_hetero = strain_hetero(aai_raw_scores, aai_strain_threshold)

    return aai_mean_bin_hetero


def strain_hetero(aai_scores, aai_strain_threshold):
    """Calculate strain heterogeneity."""
    aai_hetero = {}
    strain_count = 0
    multi_copy_pairs = 0

    for markerId in aai_scores.keys():
        local_strain_count = 0
        for aaiScore in aai_scores[markerId]:
            multi_copy_pairs += 1
            if aaiScore > aai_strain_threshold:
                strain_count += 1
                local_strain_count += 1
            hetero = float(local_strain_count) / len(aai_scores[markerId])
            aai_hetero[markerId] = hetero

    if not multi_copy_pairs == 0:
        aai_mean_bin_hetero = float(strain_count) / multi_copy_pairs
    else:
        aai_mean_bin_hetero = 0

    return aai_hetero, aai_mean_bin_hetero


def aai_seq(seq1, seq2):
    """Calculate amino acid identity between sequences."""
    assert len(seq1) == len(seq2)

    # calculation of AAI should ignore missing data at
    # the start of end of each sequence
    start_index = 0
    for i in range(0, len(seq1)):
        if seq1[i] == '-' or seq2[i] == '-':
            start_index = i + 1
        else:
            break

    end_index = len(seq1)
    for i in range(len(seq1) - 1, 0, -1):
        if seq1[i] == '-' or seq2[i] == '-':
            end_index = i
        else:
            break

    mismatches = 0
    seq_len = 0
    for i in range(start_index, end_index):
        if seq1[i] != seq2[i]:
            mismatches += 1
            seq_len += 1
        elif seq1[i] == '-' and seq2[i] == '-':
            pass
        else:
            seq_len += 1

    if seq_len == 0:
        return 0.0

    return 1.0 - (float(mismatches) / seq_len)
