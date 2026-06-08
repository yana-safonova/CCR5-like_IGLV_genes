import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from scipy import stats

def PvalToStr(pval, corr):
    if pval >= 0.05 / corr:
        return 'ns'
    if pval >= 0.01 / corr:
        return '*'
    if pval >= 0.001 / corr:
        return '**'
    if pval >= 0.0001 / corr:
        return '***'
    return '****'

species_list = ['ferret', 'sloth_bear', 'dog', 'gray_wolf', 'maned_wolf']

dir_prefix = 'data/shm_statistics'
dfs = []
for species in species_list:
    csv_path = os.path.join(dir_prefix, 'df_long_cdr1_stats_' + species + '.csv')
    df = pd.read_csv(csv_path)
    df['Species'] = species
    dfs.append(df)
df = pd.concat(dfs).reset_index(drop=True)

tissue_dict = {'PBMC' : 'PBMC', 'Spleen' : 'Spleen', 'spleen' : 'Spleen'}
df['Tissue'] = [tissue_dict[df['Tissue'][i]] for i in range(len(df))]
df['Sample'] = [df['Species'][i] + ' ' + df['Tissue'][i] for i in range(len(df))]

non_cdr1_dnds_values = []
for i in range(len(df)):
    all_n = df['All_N'][i]
    all_s = df['All_S'][i]
    cdr1_n = df['CDR1_N'][i]
    cdr1_s = df['CDR1_S'][i]
    non_cdr1_n = all_n - cdr1_n
    non_cdr1_s = all_s - cdr1_s
    non_cdr1_dnds = float('nan')
    if non_cdr1_s != 0:
        non_cdr1_dnds = non_cdr1_n / non_cdr1_s / 2
    non_cdr1_dnds_values.append(non_cdr1_dnds)
df['Dn/Ds_nonCDR1'] = non_cdr1_dnds_values

cols = ['Dn/Ds_CDR1', 'Dn/Ds_All', 'Dn/Ds_nonCDR1']
col_labels = {'CDR1' : '# SHMs in CDRL1', 'All' : '# SHMs in V gene', 'Dn/Ds_CDR1' : 'dN/dS in CDRL1', 'Dn/Ds_All' : 'dN/dS in V gene', 'Dn/Ds_nonCDR1' : 'dN/dS in non-CDRL1'}

species_colors = {'domestic_ferret': '#440154', 'black-footed_ferret': '#443983', 'red_panda': '#31688e', 'sloth_bear': '#21918c', 'domestic_dog': '#35b779', 'gray_wolf': '#90d743', 'maned_wolf': '#fde725'}
species_dict = {'ferret' : 'domestic_ferret', 'red_panda' : 'red_panda', 'sloth_bear' : 'sloth_bear', 'dog' : 'domestic_dog', 'gray_wolf' : 'gray_wolf', 'maned_wolf' : 'maned_wolf'}

palette = {}
for i in range(len(df)):
    palette[df['Sample'][i]] = species_colors[species_dict[df['Species'][i]]]

for col in cols:
    df_non_na = df.dropna(subset = [col])

    plt.figure(figsize = (8, 4))
    sns.barplot(x = 'Sample', hue = 'Sample', y = col, data = df_non_na, palette = palette)
    plt.ylabel(col_labels[col], fontsize = 14)
    plt.tight_layout()
    plt.savefig('plot_shm_stats_' + col.replace('/', '-') + '.png', dpi = 300)
    plt.clf()

    print('\n==== ' + col + ' (P-values)')
    samples = list(set(df_non_na['Sample']))
    for s1_idx, s1 in enumerate(samples):
        for s2 in samples[s1_idx + 1:]:
            h, pval = stats.kruskal(df_non_na[df_non_na['Sample'] == s1][col], df_non_na[df_non_na['Sample'] == s2][col])
            print(s1, s2, pval, PvalToStr(pval, 15))


cols = ['All', 'CDR1']
for col in cols:
    df_non_na = df[df[col] > 0]
    plt.figure(figsize = (8, 4))
    sns.barplot(x = 'Sample', hue = 'Sample', y = col, data = df_non_na, palette = palette)
    plt.ylabel(col_labels[col], fontsize = 14)
    plt.yticks(fontsize = 12)
    plt.tight_layout()
    plt.savefig('plot_shm_stats_' + col.replace('/', '-') + '.png', dpi = 300)
    plt.clf()

df_non_na = df.dropna(subset = ['Dn/Ds_CDR1', 'Dn/Ds_nonCDR1']).reset_index(drop=True)
df_non_na['DnDs_ratio'] = df_non_na['Dn/Ds_CDR1'] - df_non_na['Dn/Ds_nonCDR1']

species_colors = {'ferret': '#440154', 'black-footed_ferret': '#443983', 'red_panda': '#31688e', 'sloth_bear': '#21918c', 'dog': '#35b779', 'gray_wolf': '#90d743', 'maned_wolf': '#fde725'}
order = []
for i in range(len(df_non_na)):
    if df_non_na['Sample'][i] not in order:
        order.append(df_non_na['Sample'][i])
palette = [species_colors[o.split()[0]] for o in order]

plt.figure(figsize = (8, 4))
sns.barplot(x = 'Sample', y = 'DnDs_ratio', data = df_non_na, palette = palette)
plt.ylabel('dN/dS (CDR1) - dN/dS (non-CDR1)', fontsize = 14)
plt.xlabel('')
plt.title('long-CDRL1 IGL chains', fontsize = 14)
plt.yticks(fontsize = 12)
plt.ylim(-1.5, 8)
plt.tight_layout()
plt.savefig('plot_shm_stats_DnDs_ratio.png', dpi = 300)
plt.clf()

print('\n==== dNdS in CDR1 / non-CDRL1 (P-values)')
samples = list(set(df_non_na['Sample']))
for s1_idx, s1 in enumerate(samples):
    for s2 in samples[s1_idx + 1:]:
        h, pval = stats.kruskal(df_non_na[df_non_na['Sample'] == s1][col], df_non_na[df_non_na['Sample'] == s2][col])
        print(s1, s2, pval, PvalToStr(pval, 15))
