import os
import sys
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO

from scipy import stats
import matplotlib as mplt
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import pearsonr

import warnings
warnings.filterwarnings("ignore")

from scipy import stats

def PvalToStr(pval):
    if pval > 0.05:
        return ''
    if pval > 0.01:
        return '*'
    if pval > 0.001:
        return '**'
    if pval > 0.0001:
        return '***'
    return '****'

def ComputePValues(df, target_column, sep_column):
    states = list(set(df[target_column]))
    usage_dict = dict()
    for s in states:
        usage_dict[s] = []
    for i in range(len(df)):
        if not pd.isnull(df[sep_column][i]):
            usage_dict[df[target_column][i]].append(df[sep_column][i])
    pi_stats = ''
    if len(states) <= 1:
        return 1
    if len(states) == 2:
        pi_stats = stats.kruskal(usage_dict[states[0]], usage_dict[states[1]])
    else:
        pi_stats = stats.kruskal(usage_dict[states[0]], usage_dict[states[1]], usage_dict[states[2]])
    return pi_stats[1]

def GetColorByNormalizedValue(cmap_name, norm_value):
    if norm_value < 0 or norm_value > 1:
        print("ERROR: value " + str(norm_value) + ' does not belong to [0, 1]')
    cmap = plt.cm.get_cmap(cmap_name)
    color = cmap(norm_value)
    return mplt.colors.rgb2hex(color[:3])

def GetDistance(s1, s2, blosum_df):
    dist = 0
    for aa1, aa2 in zip(s1, s2):
        score = blosum_df.loc[blosum_df['AA'] == aa1].reset_index()[aa2][0]
        dist += score
    return np.mean(dist)

def ComputeGeneDistance(ccr5_prefix, gene, blosum_df):
    best_dist = -99
    best_pos = -1
    for i in range(len(gene) - len(ccr5_prefix) + 1):
        gene_kmer = gene[i : i + len(ccr5_prefix)]
        dist = GetDistance(gene_kmer, ccr5_prefix, blosum_df)
        if dist > best_dist:
            best_dist = dist
            best_pos = i
    return best_dist, best_pos


suborder_colors = {'Feliformia' : '#E5191D', 'Caniformia' : '#367EB8'}

family_colors = {'Canidae' : '#8ED3C7', # teal
                 'Herpestidae' : '#FB8071', # light pink
                 'Mephitidae' : '#D9D9D9', # light gray
                 'Ursidae' : '#BB80BD', # violet
                 'Eupleridae' : '#81B1D3', # blue
                 'Ailuridae' : '#BFBADA', # light violet
                 'Phocidae' : '#FCCDE5', # pink
                 'Mustelidae' : '#FDB462',  #orange
                 'Otariidae' : '#FFED6F', # yellow
                 'Felidae' : '#E5C494'} # brown


species_name_txt = 'data/species_names.txt'
species_csv = 'data/species_table_no_accessions.csv'
gene_csv = 'data/all_carnivores_IGLV.CDR1_stats-04092026.csv'
ccr5_dir = 'data/ccr5_matches'

blosum_df = pd.read_csv('data/blosum_matrix.txt')
species_names = [l.strip() for l in open(species_name_txt).readlines()]
species_df = pd.read_csv(species_csv)
gene_df = pd.read_csv(gene_csv)

out_df = {'SpeciesID' : [], 'Suborder' : [], 'Family' : [], 'CCR5Dist' : [], 'CDR1Len' : [], 'CDR1' : [], 'LatinName' : [], 'CCR5_18' : [], 'LongestCDR1' : [], 'LongestCDR1Len' : [], 'NumLongestCDR1s' : [], 'LongestCDR1s' : []}
for name in species_names:
    name_splits = name.split()
    latin_name = name #name_splits[0] + ' ' + name_splits[1]
    name_df = species_df[species_df['Species (lat)'] == name].reset_index()
    species_id = name_df['Species (eng)'][0].replace(' ', '_').lower()

    print('\n==== ' + species_id + '...')

    #### reading CCR5
    ccr5_aa_fasta = os.path.join(ccr5_dir, species_id + '_aa.fasta')
    if not os.path.exists(ccr5_aa_fasta):
        continue
    ccr5_seqs = [str(r.seq) for r in SeqIO.parse(ccr5_aa_fasta, 'fasta')]
    print('# ccr5 genes', len(ccr5_seqs))
    if len(ccr5_seqs) == 0:
        continue
    ccr5_seq = ccr5_seqs[0]
    ccr5_prefix = ccr5_seq[0 : 18]

    species_gene_df = gene_df[(gene_df['SpeciesID'] == species_id) & gene_df['Productive']].reset_index()
    max_cdr1_len = max(species_gene_df['CDR1_len'])
    max_cdr1_df = species_gene_df[species_gene_df['CDR1_len'] == max_cdr1_len].reset_index(drop = True)
    num_diff_max_cdr1s = list(set(max_cdr1_df['cdr1_aa']))
    best_dist = -100
    best_index = -1
    best_pos = -1
    for i in range(len(species_gene_df)):
        if pd.isnull(species_gene_df['cdr1_aa'][i]):
            continue
        nt_seq = species_gene_df['Sequence'][i]
        aa_seq = str(Seq(nt_seq).translate())
        cdr1_pos = aa_seq.find(species_gene_df['cdr1_aa'][i])
        dist, pos = ComputeGeneDistance(ccr5_prefix, aa_seq[cdr1_pos - 5 : cdr1_pos + species_gene_df['CDR1_len'][i] + 10], blosum_df)
        if dist > best_dist:
            best_dist = dist
            best_index = i
            best_pos = pos

    print('IGLV gene with the best CCR5 match: ' + str(Seq(species_gene_df['Sequence'][best_index]).translate()))
    print('CDRL1: ' + species_gene_df['cdr1_aa'][best_index] + ' (' + str(len(species_gene_df['cdr1_aa'][best_index])) + ' aa)')
    print('All CDR1s with max length (' + str(max_cdr1_len) + ' aa): ' + ' '.join(num_diff_max_cdr1s))

    out_df['SpeciesID'].append(species_id)
    out_df['Suborder'].append(species_gene_df['Suborder'][best_index])
    out_df['Family'].append(species_gene_df['Family'][best_index])
    out_df['LongestCDR1'].append(num_diff_max_cdr1s[0])
    out_df['LongestCDR1s'].append(','.join(num_diff_max_cdr1s))
    out_df['LongestCDR1Len'].append(max_cdr1_len)
    out_df['NumLongestCDR1s'].append(len(num_diff_max_cdr1s))
    out_df['CCR5Dist'].append(best_dist)
    out_df['CDR1Len'].append(species_gene_df['CDR1_len'][best_index])
    out_df['CDR1'].append(species_gene_df['cdr1_aa'][best_index])
    latin_name = latin_name.replace(' ', '_')
    new_latin_name = latin_name[0].upper() + latin_name[1:].lower()
    out_df['LatinName'].append(new_latin_name)
    out_df['CCR5_18'].append(ccr5_prefix)

out_df = pd.DataFrame(out_df)
out_df = out_df.sort_values(by = ['Family'])
out_df.to_csv('data/cdr1_ccr5_similarity_stats.csv')

ccr5_label = 'CDR1 / CCR5 score'

sns.barplot(x = 'Family', y = 'CCR5Dist', data = out_df, order = ['Mustelidae', 'Phocidae', 'Ursidae', 'Canidae', 'Herpestidae', 'Felidae'], palette = family_colors)
plt.xticks(rotation = 90)
plt.yticks(rotation = 90)
plt.ylabel(ccr5_label)
plt.xlabel('')
plt.tight_layout()
plt.savefig('ccr5_dist_by_family.png', dpi = 300)
plt.clf()

pval = ComputePValues(out_df, 'Suborder', 'CCR5Dist')

plt.figure(figsize = (4, 4))
sns.barplot(x = 'Suborder', y = 'CCR5Dist', data = out_df, palette = suborder_colors)
plt.xlabel('')
plt.ylabel(ccr5_label)
plt.title('P=' + PvalToStr(pval))
plt.tight_layout()
plt.savefig('ccr5_dist_by_suborder.png', dpi = 300)
plt.clf()

pval = ComputePValues(out_df, 'Suborder', 'CDR1Len')

plt.figure(figsize = (4, 4))
sns.barplot(x = 'Suborder', y = 'CDR1Len', data = out_df, palette = suborder_colors)
plt.xlabel('')
plt.ylabel('longest CDR1 length (aa)')
plt.title('P=' + PvalToStr(pval))
plt.tight_layout()
plt.savefig('cdr1_len_by_suborder.png', dpi = 300)
plt.clf()

r, pval = pearsonr(out_df['CDR1Len'], out_df['CCR5Dist'])

sns.lmplot(x = 'CDR1Len', y = 'CCR5Dist', data = out_df, height=4, line_kws = {'color' : '#7F7F7F'}, scatter_kws = {'color' : '#7F7F7F'})
plt.xlabel('CDR1 length (aa)')
plt.ylabel(ccr5_label)
plt.title('r=' + str(round(r, 2)) + ', P=' + PvalToStr(pval))
plt.tight_layout()
plt.savefig('ccr5_dist_vs_cdr1_len.png', dpi = 300)
plt.clf()

print('Correlation between CDR1 length and CCR5 similarity', r, pval)
print('Difference of CCR5 similarities between Caniformia and Feliformia', ComputePValues(out_df, 'Suborder', 'CCR5Dist'))
