import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np
import itertools

def PvalToStr(pval):
    if pval >= 0.05:
        return 'ns'
    if pval >= 0.01:
        return '*'
    if pval >= 0.001:
        return '**'
    if pval >= 0.0001:
        return '***'
    return '****'

input_dir = 'data/repertoire_gene_usages'
csv_files = [f for f in os.listdir(input_dir) if f.find('df_cdr1') != -1]
dfs = [pd.read_csv(os.path.join(input_dir, csv)) for csv in csv_files]
df = pd.concat(dfs).reset_index()

species = set(df['Species'])

species_dict = {'red_panda' : 'Red panda', 'sloth_bear' : 'Sloth bear', 'ferret' : 'Domestic ferret', 'maned_wolf' : 'Maned wolf', 'bff' : 'Black-footed ferret', 'dog' : 'Domestic dog', 'gray_wolf' : 'Gray wolf'}
species_order = ['ferret', 'bff', 'red_panda', 'sloth_bear', 'dog', 'gray_wolf', 'maned_wolf']
tissue_dict = {'PBMC' : 'PBMC', 'spleen' : 'Spleen', 'Spleen' : 'Spleen'}

df['Tissue'] = [tissue_dict[df['Tissue'][i]] for i in range(len(df))]

cdr1_lens = sorted(set(df['Len']))

long_cdr1_df = {'Species' : [], 'Tissue' : [], 'LongCDR1' : [], 'Sample' : []}
subjects = set(df['Subject'])
for subj in subjects:
    subj_df = df[df['Subject'] == subj]
    species = subj_df['Species'].iloc[0]
    subj_tissues = set(subj_df['Tissue'])
    for tissue in subj_tissues:
        long_sub_df = subj_df[(subj_df['Len'] >= 13) & (subj_df['Tissue'] == tissue)].reset_index()
        long_cdr1_df['Species'].append(species)
        long_cdr1_df['Tissue'].append(tissue)
        long_cdr1_df['Sample'].append(species + ' ' + tissue)
        perc_long_cdr1 = 0
        if len(long_sub_df) != 0:
            perc_long_cdr1 = sum(long_sub_df['Percentage'])
        long_cdr1_df['LongCDR1'].append(perc_long_cdr1)

order = []
for species in species_order:
    for t in sorted(set(df['Tissue'])):
        sub_df = df[(df['Species'] == species) & (df['Tissue'] == t)]
        if len(sub_df) != 0:
            order.append(species + ' ' + t)

species_colors = {'ferret': '#440154', 'bff': '#443983', 'red_panda': '#31688e', 'sloth_bear': '#21918c', 'dog': '#35b779', 'gray_wolf': '#90d743', 'maned_wolf': '#fde725'}
palette = [species_colors[o.split()[0]] for o in order]

long_cdr1_df = pd.DataFrame(long_cdr1_df)
sub_df = long_cdr1_df[long_cdr1_df['Species'].isin(['ferret', 'dog'])]

plt.figure(figsize = (12, 3))
sns.barplot(x = 'Sample', y = 'LongCDR1', data = long_cdr1_df, order = order, palette = palette, err_kws={'linewidth': 1})
sns.swarmplot(x = 'Sample', y = 'LongCDR1', data = sub_df, order = order, color = 'black', edgecolor = 'white', linewidth = 0.5)
plt.xlabel('')
plt.ylabel('% long-CDRL1 IGL chains')
plt.tight_layout()
plt.savefig('plot_long_cdr1_fractions.png', dpi = 300)
plt.clf()

print('Avg usages of long-CDRL1 IGLV genes')
for species in set(long_cdr1_df['Species']):
    s_df = long_cdr1_df[long_cdr1_df['Species'] == species]
    print(' ', species, np.mean(s_df['LongCDR1']))

group_dict = {'ferret' : 'Arctoidea', 'bff' : 'Arctoidea', 'red_panda' : 'Arctoidea', 'sloth_bear' : 'Arctoidea', 'dog' : 'Canidae', 'gray_wolf' : 'Canidae', 'maned_wolf' : 'Canidae'}
long_cdr1_df['Group'] = [group_dict[long_cdr1_df['Species'][i]] for i in range(len(long_cdr1_df))]
long_cdr1_df = long_cdr1_df[long_cdr1_df['Species'] != 'bff']
long_cdr1_df = long_cdr1_df.groupby(['Species']).agg({'Species' : 'first', 'Group' : 'first', 'LongCDR1' : 'mean'})

t, pval = stats.kruskal(long_cdr1_df[long_cdr1_df['Group'] == 'Arctoidea']['LongCDR1'], long_cdr1_df[long_cdr1_df['Group'] == 'Canidae']['LongCDR1'])
print('\nUsage of long-CDRL1 genes in Arctoidea vs Canidae', pval)
palette = [species_colors[long_cdr1_df['Species'][i]] for i in range(len(long_cdr1_df))]

plt.figure(figsize = (4, 4))
sns.barplot(x = 'Group', y = 'LongCDR1', data = long_cdr1_df, color = '#C7C7C7', err_kws={'linewidth': 1}, order = ['Arctoidea', 'Canidae'])
ax = sns.swarmplot(x = 'Group', y = 'LongCDR1', hue = 'Species', data = long_cdr1_df, palette = palette, order = ['Arctoidea', 'Canidae'])
for col in ax.collections:
    col.set_sizes([75])
ax.get_legend().remove()
plt.xlabel(' ')
plt.ylabel('% long-CDRL1 IGL chains', fontsize = 14)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.title('P=' + PvalToStr(pval), fontsize = 14)
plt.tight_layout()
plt.savefig('plot_Arctoidea_vs_Canidae.png', dpi = 300)
plt.clf()
