# Comparative study of CCR5-like IGLV genes across carnivorans

## Annotation of IGLV and CCR5 across carnivorans
- Annotations of productive IGLV genes (in-frame IGLV genes without stop codons) across 47 carnivorans are combined in `data/all_carnivores_IGLV.CDR1_stats-04092026.csv`. It includes sequences of IGLV genes, their positions in the genome assemblies, and the corresponding CDRL1 sequences.

- Annotations of CCR5 gene matches are collected in `data/ccr5_matches`.

## Similarity of antigen-binding sites of IGLV and CCR5 genes
- To compute similarity between CCR5 and IGLV binding sites for the domestic ferret, run:
  
`python scripts/compute_ccr5_similarity_ferret.py`

- To compute similarity between CCR5 and IGLV binding sites across 47 carnivorans, run:

`python scripts/compute_cdrl1_ccr5_similarity.py`

The script reports the dataframe `data/cdr1_ccr5_similarity_stats.csv` containing similarity statistics and plots visualizing them: `ccr5_dist_by_family.png`, `ccr5_dist_by_suborder.png`, `ccr5_dist_vs_cdr1_len.png`, `cdr1_len_by_suborder.png`.

## dN/dS statistics across caniforms 
- dN/dS statistics for species across four Caniformia families (Canidae, Mustelidae, Phocidae, Ursidae) are collected in `data/germline_dnds`. Each file contains the following columns:

  - `target_sample`, `target_gene`, `target_aa`, `target_cdr1` - a IGLV gene in species `target_sample` with its nucleotide sequence, amino acid sequence, and CDRL1.  
  - `matched_sample`, `matched_gene`, `matched_aa` - the closest IGLV gene in species `matched_sample`.
  - `IsTargetGene` - whether the target IGLV gene is long-CDRL1.
  - `dist` - percent identity between `target_gene` and `matched_gene`.
  - `NumS`, `NumN`, `dNdS` - counts of synonymous and non-synonymous mutations as well as the corresponding dN/dS value. 

- To compute and visualize the differences between dN/dS values for long-CDRL1 and short-CDRL1 IGLV genes across four Caniformia families, run:

`python scripts/visualize_germline_dnds.py`

## Antibody repertoire characteristics 
- Antibody repertoire usages of IGLV genes broken down by CDRL1 lengths of V(D)J recombinations across seven carnivorans are collected in `data/repertoire_gene_usages`.

- To visualize usages of long-CDRL1 IGLV genes, run:

`python scripts/visualize_repertoire_usages.py`

The script also computes the difference in usages between Arctoidea and Canidae species.

- SHM characteristics of long-CRL1 IGLV genes are collected in `data/shm_statistics`. Each file includes the following columns:
  
  - `All`, `All_N`, `All_S` - counts of SHMs (all, non-synonymous only, synonymous only) in the region spanning FR1-FR3.
  - `CDR1`, `CDR1_N`, `CDR1_S` -  - counts of SHMs (all, non-synonymous only, synonymous only) in CDRL1.

- To visualize SHM characteristics in `data/shm_statistics`, run:
  
`python scripts/visualize_shm_stats.py`

## Citation
Safonova Y, Pursell T, Sheneman KR, Whitley CS, Mikhailova A, Pattar V, Pospelova M, Rubio AA, Voss KA, Welker JM, Zamyatin A, Bankevich A, Boeke JD, Haraguchi E, Hudson E, Kline E, Lama TM, Lauer W, Le Sage V, Thomas M, Watson CT, Zheng S, Barnes CO, Lakdawala SS, Pennell M, Smith ML, Boyd S, Lawrenz MB, Koepfli K-P. Disruption of a CCR5-like immunoglobulin gene is linked to plague susceptibility in black-footed ferrets. To be submitted. 2026.
