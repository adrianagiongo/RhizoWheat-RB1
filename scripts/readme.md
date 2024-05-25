### Github repository for 

## Dataset jki_seq5 (16S rRNA gene amplicon sequencing) pipeline analyses
Performed in R v.4.1.3 [R Core Team](https://www.r-project.org)

#### Color code
- Depth
  - ![#70411c](https://placehold.co/15x15/70411c/70411c.png) `#70411c` (D30)
  - ![#ba997b](https://placehold.co/15x15/ba997b/ba997b.png) `#ba997b` (D60)
- Microhabitats
  - ![#d9b967](https://placehold.co/15x15/d9b967/d9b967.png) `#d9b967` (RH)
  - ![#57896A](https://placehold.co/15x15/57896A/57896A.png) `#57896A` (RP)
- Rotations
  - ![#d7d3a9](https://placehold.co/15x15/d7d3a9/d7d3a9.png) `#d7d3a9` (W1)
  - ![#74a553](https://placehold.co/15x15/74a553/74a553.png) `#74a553` (W2)

### Packages
library("dada2")\
library("ShortRead")\
library("Biostrings")\
library("phyloseq")\
library("microbiome")\
library("vegan")\
library("DESeq2") \
library("dplyr")\
library("stringr")\
library("ggpubr")\
library("tidyr")\
library("ggplot2")\
library("readxl")\
library("RColorBrewer")
  
### 1. Dada2
This script uses the raw data obtained from the BioProject SRA.\
A R server is required.\
Database used: [SILVA 138 SSU](https://www.arb-silva.de/documentation/release-138/)
  - Non-truncate parameter
```
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, compress = TRUE, truncQ = 2, trimLeft=c(19,20), minLen = 100, maxN=0, maxEE=c(2,2), rm.phix=TRUE, matchIDs = TRUE, multithread=TRUE)
```
Total of samples = 30\
Total of reads (nonchim) = 2547216

### 2. Creating phyloseq object
This script creates a phyloseq object (package phyloseq) based on these files:
- jki_seq1_metadata.csv
- jki_seq1_otu.xlsx
- jki_seq1_seqs.xlsx
- jki_seq1_taxa.xlsx

Initial number of samples = 36 (18284 ASVs)

### 3. Rename NA
This script replaces NA for the latest taxonomy found for an ASV.\
Source = https://github.com/joey711/phyloseq/issues/850

### 4. Clean dataset
Initial number of samples = 36 (18284 ASVs)

This script removes unwanted taxonomic groups from the dataset.
- Root NAs (Domain, Phylum)
- Eukaryotes (Domain) (42 ASVs removed)
- Chloroplasts (Order) (47 ASVs removed)
- Mitochondria (Family) (21 ASVs removed)

Total of samples = 36 (18174 ASVs)

### 5. Rarefaction
This script performs rarefaction based on the minimum sequences. See Schloss (2024).

2312 ASVs removed after rarefaction.\
Final = 15862 ASVs and 34391 sequences

### 6. Data selection
Based on the rarefied data, this script group taxa by desired rank, selects a group of samples to be analyzed separately and filter them to remover zeros.

Phyla = 47\
Annotation = 1326\
Besides 'all', samples were grouped by T1 and T2

T1 = 11960 ASVs / 24 samples\
T2 = 7973 ASVs / 12 samples

### 7. Alpha diversity
This script calculates alpha diversity based on the rarefied data.

No significantly differences in Shannon, Chao1 and Pielou (Wilcoxon test)

### 8. Ordination 
This script creates MDS plots and calculates PERMANOVA and ANOSIM based on the rarefied data. Square root was applied (*hellinger*).\

#### PERMANOVA

| PERMANOVA All               | Df | SumOfSqs | R2     | F      | Pr(>F) |
| --------------------------- | -- | -------- | ------ | ------ | ------ |
| Rotation                    | 1  | 0.1855   | 0.0531 | 2.242  | 0.0001 |
| Depth                       | 1  | 0.2735   | 0.0784 | 3.3057 | 0.0001 |
| Microhabitat                | 1  | 0.1612   | 0.0462 | 1.9482 | 0.0004 |
| Stage                       | 1  | 0.236    | 0.0676 | 2.8529 | 0.0001 |
| Rotation:Depth              | 1  | 0.1039   | 0.0298 | 1.2558 | 0.0765 |
| Rotation:Microhabitat       | 1  | 0.0854   | 0.0245 | 1.0323 | 0.3515 |
| Depth:Microhabitat          | 1  | 0.0967   | 0.0277 | 1.1682 | 0.142  |
| Rotation:Stage              | 1  | 0.1064   | 0.0305 | 1.2863 | 0.0603 |
| Depth:Stage                 | 1  | 0.1044   | 0.0299 | 1.2613 | 0.0709 |
| Rotation:Depth:Microhabitat | 1  | 0.0741   | 0.0212 | 0.895  | 0.7196 |
| Rotation:Depth:Stage        | 1  | 0.078    | 0.0224 | 0.9432 | 0.582  |
| Residual                    | 24 | 1.9858   | 0.5688 | NA     | NA     |
| Total                       | 35 | 3.4909   | 1      | NA     | NA     |

| PERMANOVA T1                | Df | SumOfSqs | R2     | F      | Pr(>F) |
| --------------------------- | -- | -------- | ------ | ------ | ------ |
| Microhabitat                | 1  | 0.1172   | 0.0523 | 1.3327 | 0.0344 |
| Rotation                    | 1  | 0.1593   | 0.0712 | 1.8122 | 0.0006 |
| Depth                       | 1  | 0.2223   | 0.0993 | 2.5284 | 0.0001 |
| Microhabitat:Rotation       | 1  | 0.0711   | 0.0317 | 0.8084 | 0.9252 |
| Microhabitat:Depth          | 1  | 0.0915   | 0.0409 | 1.0410 | 0.3433 |
| Rotation:Depth              | 1  | 0.1002   | 0.0447 | 1.1391 | 0.1666 |
| Microhabitat:Rotation:Depth | 1  | 0.0708   | 0.0316 | 0.8048 | 0.9308 |
| Residual                    | 16 | 1.4067   | 0.6283 | NA     | NA     |
| Total                       | 23 | 2.2390   | 1.0000 | NA     | NA     |

| PERMANOVA T2   | Df | SumOfSqs | R2     | F      | Pr(>F) |
| -------------- | -- | -------- | -------| ------ | -------|
| Rotation       | 1  | 0.1469   | 0.1512 | 2.0302 | 0.0005 |
| Depth          | 1  | 0.1607   | 0.1654 | 2.2203 | 0.0001 |
| Rotation:Depth | 1  | 0.0851   | 0.0875 | 1.1755 | 0.1751 |
| Residual       | 8  | 0.5790   | 0.5959 | NA     | NA     |
| Total          | 11 | 0.9718   | 1      | NA     | NA     |

#### anosim
Permutation: free
Number of permutations: 10000

All: ANOSIM statistic R: 0.3859 
- Significance: 9.999e-05 

T1: ANOSIM statistic R: 0.1135 
- Significance: 0.030997 

T2: ANOSIM statistic R: 0.4426 
- Significance: 0.0083992 

### 9. Heatmap
Heatmaps based on the Top10 taxa (Annotation) based on the rarefied data. Relative abundance or square root (hellinger).

### 10. Pheatmap
Pheatmaps of the Top10 taxa (Annotation) based on the rarefied data. Relative abundance (hellinger).\
Column clustered by Depth, Microhabitat and Rotation. T1 and T2.

### 9. DESeq2
Based on the rarefied data, this script performs differential abundance (DA) between two groups.

### References
Callahan B, McMurdie P, Rosen M et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat. Met. 13:581–583. [DOI](http://10.1101/024034)

Schloss PD (2024) Rarefaction is currently the best approach to control for uneven sequencing effort in amplicon sequence analyses. mSphere 9:e00354-23. [DOI](http://10.1128/msphere.00354-23)

