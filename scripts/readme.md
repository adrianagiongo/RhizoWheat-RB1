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





| PERMANOVA All   | Df  | SumOfSqs | R2      | F      | Pr(>F)  |
| --------------- |:---:| --------:| -------:| ------ | ------- |
| Rotation        | 3   | 0.8443   | 0.15947 | 1.8312 | 0.03110 |
| Depth           | 1   | 0.4626   | 0.08739 | 3.0105 | 0.01550 |
| Stage           | 3   | 0.7600   | 0.14355 | 1.6484 | 0.05509 |
| Microhabitat    | 3   | 0.7600   | 0.14355 | 1.6484 | 0.05509 |
| Microhabitat    | 3   | 0.7600   | 0.14355 | 1.6484 | 0.05509 |
| Microhabitat    | 3   | 0.7600   | 0.14355 | 1.6484 | 0.05509 |
| Residual        | 21  | 3.2272   | 0.60959 |        |         |
| Total           | 29  | 5.2941   | 1.00000 |        |         |

| PERMANOVA T1    | Df  | SumOfSqs | R2      | F      | Pr(>F)  |
| --------------- |:---:| --------:| -------:| ------ | ------- |
| Treatment       | 3   | 0.30482  | 0.24697 | 0.9839 | 0.6315  |
| Residual        | 9   | 0.92942  | 0.75303 |        |         |
| Total           | 12  | 1.23424  | 1.00000 |        |         |

| PERMANOVA T2    | Df  | SumOfSqs | R2      | F      | Pr(>F)  |
| --------------- |:---:| --------:| -------:| ------ | ------- |
| Treatment       | 3   | 1.3195   | 0.36478 | 2.297  | 0.0126  |
| Residual        | 12  | 2.2978   | 0.63522 |        |         |
| Total           | 15  | 3.6173   | 1.00000 |        |         |

### 9. DESeq2
Based on the rarefied data, this script performs differential abundance (DA) between two groups.

### References
Callahan B, McMurdie P, Rosen M et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat. Met. 13:581–583. [DOI](http://10.1101/024034)

Schloss PD (2024) Rarefaction is currently the best approach to control for uneven sequencing effort in amplicon sequence analyses. mSphere 9:e00354-23. [DOI](http://10.1128/msphere.00354-23)

