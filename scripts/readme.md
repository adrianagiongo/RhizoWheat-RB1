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
  - ![#9cc184](https://placehold.co/15x15/9cc184/9cc184.png) `#9cc184` (W1)
  - ![#447243](https://placehold.co/15x15/447243/447243.png) `#447243` (W2)

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
A R server is required. \
Database used: [SILVA 138 SSU](https://www.arb-silva.de/documentation/release-138/) 
- Non-truncate parameter

### 2. Creating phyloseq object
This script creates a phyloseq object based on these files:

- jki_seq1_metadata.csv
- jki_seq1_otu.xlsx
- jki_seq1_seqs.xlsx
- jki_seq1_taxa.xlsx

### 3. Rename NA
This script replaces NA for the latest taxonomy found for an ASV.

### 4. Clean dataset
This script removes unwanted taxonomic groups from the dataset.
- Root NAs (Domain, Phylum)
- Eukaryotes (Domain)
- Chloroplasts (Order)
- Mitochondria (Family)

### 5. Data selection
This script selects a group of samples to be analyzed separately.

### 6. Rarefaction
This script performs rarefaction based on the minimum sequences.

### 7. Alpha diversity
This script calculates alpha diversity based on the rarefied data.

### 8. Ordination 
This script creates MDS plots and calculates PERMANOVA and ANOSIM based on the rarefied data.

### 9. DESeq2
Based on the rarefied data, this script performs differential abundance (DA) between two groups.
