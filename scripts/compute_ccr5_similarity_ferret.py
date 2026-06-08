import os
import sys
import pandas as pd
import random
import numpy as np
from scipy import stats

import matplotlib.pyplot as plt
import seaborn as sns

from Bio import SeqIO
from Bio.Seq import Seq

def GetScore(aa1, aa2, df):
    if aa1 == '-' or aa2 == '-':
        return -4
    return df[df['AA'] == aa1].reset_index()[aa2].iloc[0]

def ComputeScore(seq1, seq2, df):
    return sum([GetScore(aa1, aa2, df) for aa1, aa2 in zip(seq1, seq2)])

def ComputeKmerDistances(iglv_gene, ccr5_18, df):
    scores = []
    for i in range(len(iglv_gene) - len(ccr5_18) + 1):
        score = ComputeScore(iglv_gene[i : i + len(ccr5_18)], ccr5_18, df)
        scores.append(score)
    return scores

def GenerateRandomAA(seq_len):
    aa_list = 'A R N D C Q E G H I L K M F P S T W Y V'.split()
    aa_seq = []
    for i in range(seq_len):
        aa_seq.append(aa_list[random.randint(0, len(aa_list) - 1)])
    return ''.join(aa_seq)

def ComputeRandomScores(ccr5_18, df, num_iter = 1000):
    seq_len = len(ccr5_18)
    rand_scores = []
    for i in range(num_iter):
        rand_seq = GenerateRandomAA(seq_len)
        rand_scores.append(ComputeScore(ccr5_18, rand_seq, df))
    return rand_scores


iglv_gene = 'WAVLTQPPYVSGALGESVTISCTGIPTSIDYDEEEYTYNVNWYQQLQGKVPILLIYGDNNRNPGVPDRFSGSKSGSSASLTISGLQAEDEADYYCQSTDSSFSA'
cdr1_18 = 'TGIPTSIDYDEEEYTYNV'
ccr5_18 = 'MNDPTPTPYYDIDYGMSE'
ccr5_prefix = '                      MNDPTPTPYYDIDYGMSEPCQKTDVRKIAARLLPPLYSLVFVFGFVGNMLVVLILINCKRLKSMTDIYLLNLAISDLLL'.replace(' ', '-')
human_gene = 'QSVLTQPPSVSGAPGQRVTISCTGSSSNIGAG-----YDVHWYQQLPGTAPKLLIYGNSNRPSGVPDRFSGSKSGTSASLAITGLQAEDEADYYCQSYDSSLSG'

blosum_txt = 'data/blosum_matrix.txt'
blosum_df = pd.read_csv(blosum_txt)

iglv_genes = [r for r in SeqIO.parse('data/germline_genes/domestic_ferret_IGLV.fa', 'fasta')]
aa_seqs = []
for r in iglv_genes:
    if r.id == 'IGLV1-73*01':
        continue
    aa_seqs.append(str(r.seq.translate()))

other_iglv_scores = []
for aa_seq in aa_seqs:
    cur_scores = ComputeKmerDistances(aa_seq, ccr5_18, blosum_df)
    other_iglv_scores.extend(cur_scores)
other_iglv_mean_score = np.mean(other_iglv_scores)

kmer_scores = ComputeKmerDistances(iglv_gene, ccr5_18, blosum_df)
sorted_scores = sorted(kmer_scores)

best_idx = iglv_gene.find(cdr1_18)
print('CCR5 binding site: ' + ccr5_18)
print('Best long-CDRL1 IGLV match: ' + iglv_gene[best_idx : best_idx + len(ccr5_18)])

random_scores = ComputeRandomScores(ccr5_18, blosum_df, 5000)
mean_random_score = np.mean(random_scores)

best_score = kmer_scores[best_idx]
other_scores = kmer_scores[ : best_idx] + kmer_scores[best_idx + 1 :]
mean_other_score = np.mean(other_scores)

print('\n== P-values for the best match:')
print('vs long-CDRL1 k-mers',  np.mean(np.abs(other_scores - mean_other_score) >= abs(best_score - mean_other_score)))
print('vs random k-mers', np.mean(np.abs(random_scores - mean_random_score) >= abs(best_score - mean_random_score)))
print('vs other ferret IGLV k-mers', np.mean(np.abs(other_iglv_scores - other_iglv_mean_score) >= abs(best_score - other_iglv_mean_score)))

print('\nP-value IGLV k-mers vs random', stats.kruskal(other_scores, random_scores)[1])

print('\n== Scores:')
print('long-CDRL1 score', best_score)
print('other long-CDRL1 18-mers score', mean_other_score)
print('other IGLV 18-mer score', other_iglv_mean_score)
print('random 18-mer score', mean_random_score)

sum_df = {'Score' : other_scores + [best_score] + random_scores + other_iglv_scores, 'Type' : ['other 18-mers'] * (len(sorted_scores) - 1) + ['CDR1 18-mer'] + ['random 18-mers'] * len(random_scores) + ['other IGLV 18-mers'] * len(other_iglv_scores)}
sum_df = pd.DataFrame(sum_df)

plt.figure(figsize = (6, 4))
sns.barplot(x = 'Type', y = 'Score', data = sum_df, order = ['CDR1 18-mer', 'other 18-mers', 'other IGLV 18-mers', 'random 18-mers'], palette = ['#9567BD', '#AD8CC9', '#C4B0D4', '#C7C7C7'])
plt.xlabel('')
plt.ylabel('total BLOSUM score', fontsize = 14)
plt.xticks([0, 1, 2, 3], ['long-CDRL1', 'other 18-mers\nlong-CDRL1 IGLV', 'all 18-mers\nother IGLVs', 'random'], fontsize = 10)
plt.ylim(-20, 20)
plt.yticks(fontsize = 12)
plt.tight_layout()
plt.savefig('plot.png', dpi = 300, transparent = True)
