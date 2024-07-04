### GitHub repository for 

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
Total of samples = 36\
Total of reads (nonchim) = 1,668,656

|     | input  | filtered | denoisedF | denoisedR | merged | nonchim |
| --- | ------ | -------- | --------- | --------- | ------ | ------- |
| D01 | 86392  | 73734    | 68819     | 67865     | 43371  | 36420   |
| D02 | 93073  | 79765    | 75073     | 74316     | 49816  | 42091   |
| D03 | 95502  | 86134    | 81688     | 81091     | 59381  | 53307   |
| D04 | 96681  | 87727    | 82343     | 81377     | 54852  | 48174   |
| D05 | 80813  | 73098    | 68004     | 67241     | 40531  | 34414   |
| D06 | 95956  | 87204    | 80621     | 79865     | 47557  | 39520   |
| D07 | 80723  | 73087    | 68276     | 67310     | 42208  | 36043   |
| D08 | 119603 | 108426   | 101165    | 100306    | 61147  | 50698   |
| D09 | 97046  | 87890    | 82216     | 81508     | 52561  | 44395   |
| D10 | 110285 | 99725    | 92439     | 91871     | 59829  | 51004   |
| D11 | 112856 | 102394   | 96269     | 95293     | 63512  | 54322   |
| D12 | 88439  | 80433    | 75079     | 73838     | 46057  | 39584   |
| D13 | 83557  | 75589    | 70656     | 69996     | 44678  | 37811   |
| D14 | 112411 | 101890   | 95906     | 95319     | 63668  | 54514   |
| D15 | 110759 | 100141   | 94807     | 94299     | 68088  | 60853   |
| D16 | 102340 | 93440    | 88028     | 86769     | 59088  | 51716   |
| D17 | 91119  | 82361    | 76781     | 75762     | 46219  | 38743   |
| D18 | 96362  | 87447    | 81999     | 80988     | 52041  | 44080   |
| D19 | 86523  | 78316    | 73226     | 72603     | 46268  | 39603   |
| D20 | 93012  | 84086    | 78059     | 77223     | 45618  | 38571   |
| D21 | 86401  | 78114    | 72984     | 72310     | 46756  | 40280   |
| D22 | 102721 | 93032    | 86911     | 86499     | 54701  | 46991   |
| D23 | 100767 | 91228    | 85608     | 84866     | 55921  | 48760   |
| D24 | 85409  | 77435    | 72371     | 71674     | 47686  | 42238   |
| D25 | 108026 | 98057    | 92040     | 91165     | 57730  | 48660   |
| D26 | 107024 | 97180    | 91373     | 90914     | 60737  | 53231   |
| D27 | 95005  | 86215    | 81439     | 80766     | 55791  | 49705   |
| D28 | 113361 | 103394   | 97392     | 96830     | 65158  | 57346   |
| D29 | 109020 | 98334    | 91431     | 90061     | 51767  | 42523   |
| D30 | 103886 | 94399    | 88146     | 87064     | 54342  | 46269   |
| D31 | 85111  | 76971    | 71679     | 70855     | 45134  | 39036   |
| D32 | 110473 | 99530    | 92356     | 91537     | 54806  | 46958   |
| D33 | 103093 | 93249    | 87269     | 86332     | 55233  | 47722   |
| D34 | 111829 | 101017   | 94547     | 93520     | 60334  | 52823   |
| D35 | 111005 | 100288   | 93962     | 93287     | 61132  | 53160   |
| D36 | 118231 | 107276   | 100693    | 99368     | 65098  | 57091   |

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

### 5. Rarefaction
This script performs rarefaction based on the minimum sequences. See Schloss (2024).

2312 ASVs removed after rarefaction.\
Final = 15862 ASVs and 34391 sequences

### 6. Data selection
Based on the rarefied data, this script group taxa by desired rank, selects a group of samples to be analyzed separately and filter them to remover zeros.

All
  - Phyla = 47 / Annotations = 1326\

Samples were grouped by:
Stage (T1 and T2)
  - T1 = 11960 ASVs / 24 samples / 1202 Annotations
  - T2 = 7973 ASVs / 12 samples / 980 Annotations

Stage (T1 and T2) and Depth (D30 and D60)
  - T1_D30 = 7738 ASVs / 12 samples / 949 Annotations
  - T2_D30 = 4973 ASVs / 6 samples / 755 Annotations
  - T1_D60 = 7981 ASVs / 12 samples / 999 Annotations
  - T2_D60 = 5525 ASVs / 6 samples / 807 Annotations
    
### 7. Alpha diversity
This script calculates alpha diversity based on the rarefied data.

No significant differences in Shannon, Chao1, and Pielou (Wilcoxon test) in W1 vs W2 in T1 or T2.


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
Heatmaps based on the Top10 taxa (Annotation) based on the rarefied data. Relative abundance or square root (*hellinger*).

### 10. Pheatmap
Pheatmaps of the Top10 taxa (Annotation) based on the rarefied data. Relative abundance (*hellinger*).\
Column clustered by Depth, Microhabitat and Rotation. T1 and T2.

### 11. DESeq2
Based on the rarefied data, this script performs differential abundance (DA) between two groups.

Comparisons (*p* adjusted)
  - T1
      - W1 vs W2 (16 DA taxa) &#8594; (Bacteroidota, Firmicutes, Patescibacteria, Proteobacteria, Verrucomicrobiota)
  - T2
      - W1 vs W2 (17 DA taxa) &#8594; (Bacteroidota, Firmicutes, Patescibacteria, Proteobacteria, Verrucomicrobiota)
   
  - T1_D30
      - W1 vs W2 (12 DA taxa) &#8594; (Actinobacteriota, Firmicutes, Patescibacteria, Proteobacteria)
  - T2_D30
      - W1 vs W2 (5 DA taxa) &#8594; (Acidobacteriota, Proteobacteria)
   
  - T1_D60
      - W1 vs W2 (16 DA taxa) &#8594; (Acidobacteriota, Bacteroidota, Chloroflexi, Firmicutes, Cyanobacteria, Myxococcota, Patescibacteria, Proteobacteria, Verrucomicrobiota)
  - T2_D60
      - W1 vs W2 (5 DA taxa) &#8594; (Acidobacteriota, Proteobacteria)
  
### 12. Bubble plot (DESeq2 results)
Ggplot2 based on the DESeq2 results.

### 13. ASVs to FASTA (DESeq2 results)
Based on the DESeq2 results, we selected the ASVs from each group (T1 D30, T1 D60, T2 D30, T2 D60) using this script. After, we ran the fasta containing the ASVs against whole genomes.

### 14. ASVs to FASTA (for Tax4Fun2)

### 15. Tax4Fun2

### References
Callahan B, McMurdie P, Rosen M et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat. Met. 13:581–583. [DOI](http://10.1101/024034)

Schloss PD (2024) Rarefaction is currently the best approach to control for uneven sequencing effort in amplicon sequence analyses. mSphere 9:e00354-23. [DOI](http://10.1128/msphere.00354-23)

