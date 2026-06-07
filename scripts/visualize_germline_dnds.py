import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np

def p_value_to_stars(p_value):
    if p_value > 0.05:
        return "ns"
    if p_value <= 0.0001:
        return "****"
    if p_value <= 0.001:
        return "***"
    if p_value <= 0.01:
        return "**"
    return "*"

data_dir = 'data/germline_dnds'
csv_files = ['summary_Canidae.csv', 'summary_Mustelidae.csv', 'summary_Phocidae.csv', 'summary_Ursidae.csv']
families = [csv.split('.')[0].split('_')[1] for csv in csv_files]

hue_palette = {'long-CDRL1' : '#825CA6', 'short-CDRL1' : '#C5B0D6'}
hue_dict = {True : 'long-CDRL1', False : 'short-CDRL1'}

for csv, family in zip(csv_files, families):
    df = pd.read_csv(os.path.join(data_dir, csv))
    df['IsTargetGene'] = [hue_dict[df['IsTargetGene'][i]] for i in range(len(df))]
    df = df.dropna(subset = ['dNdS'])
    r, pval = stats.kruskal(df[df['IsTargetGene'] == 'long-CDRL1']['dNdS'], df[df['IsTargetGene'] != 'long-CDRL1']['dNdS'])

    print('====')
    print(family + ', P = ' +  str(pval))
    print('avg dN/dS (long-CDRL1 IGLV): ' + str(np.mean(df[df['IsTargetGene'] == 'long-CDRL1']['dNdS'])) + ', avg dN/dS (short-CDRL1 IGLV): ' + str(np.mean(df[df['IsTargetGene'] != 'long-CDRL1']['dNdS'])))
    plt.figure(figsize = (4, 4))
    sns.barplot(x = 'IsTargetGene', hue = 'IsTargetGene', y = 'dNdS', data = df, order = ['long-CDRL1', 'short-CDRL1'], palette = hue_palette)
    plt.xlabel(' ')
    plt.ylabel('dN/dS', fontsize = 14)
    plt.title(family + ', P=' + p_value_to_stars(pval), fontsize = 14)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    plt.ylim(0, 3.2)
    plt.tight_layout()
    plt.savefig('plot_dnds_' + family + '.png', dpi = 300)
    plt.clf()
    
