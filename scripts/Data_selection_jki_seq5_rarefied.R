###Select groups and subset data from the cleaned dataset
##use on a phyloseq object (large phyloseq data)
#load packages
library("phyloseq")
library("dplyr")
library("microbiome")

#Calculate number of genera present after filtration
length(get_taxa_unique(psO_jki_seq5_rarefied, taxonomic.rank = "Phylum"))
length(get_taxa_unique(psO_jki_seq5_rarefied, taxonomic.rank = "Annotation"))

#Agglomerate taxa using tax glom function and removing NAs
psO_jki_seq5_rarefied_phylum <- tax_glom(psO_jki_seq5_rarefied, "Phylum", NArm= TRUE)
psO_jki_seq5_rarefied_phylum
psO_jki_seq5_rarefied_annotation <- tax_glom(psO_jki_seq5_rarefied, "Annotation", NArm= TRUE)
psO_jki_seq5_rarefied_annotation
ntaxa(psO_jki_seq5_rarefied); ntaxa(psO_jki_seq5_rarefied_annotation)

#### Create tables 
df_psO_jki_seq5_rarefied_phylum <- data.frame(tax_table(psO_jki_seq5_rarefied_phylum),otu_table(psO_jki_seq5_rarefied_phylum))
write.csv(df_psO_jki_seq5_rarefied_phylum, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/df_psO_jki_seq5_rarefied_phylum.csv")

df_psO_jki_seq5_rarefied_annotation <- data.frame(tax_table(psO_jki_seq5_rarefied_annotation),otu_table(psO_jki_seq5_rarefied_annotation))
write.csv(df_psO_jki_seq5_rarefied_annotation, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/df_psO_jki_seq5_rarefied_annotation.csv")

###Transforming absolute abundance to relative abundance
##Using transform_sample_counts()
#relative abundance (%)

#For psO_jki_seq5_rarefied_phylum
psO_jki_seq5_rarefied_phylum_rel<-transform_sample_counts(psO_jki_seq5_rarefied_phylum, function(x) (x*100)/sum(x))
psO_jki_seq5_rarefied_phylum_rel
head(otu_table(psO_jki_seq5_rarefied_phylum_rel))

#### Create tables 
df_psO_jki_seq5_rarefied_phylum_rel <- data.frame(tax_table(psO_jki_seq5_rarefied_phylum_rel),otu_table(psO_jki_seq5_rarefied_phylum_rel))
write.csv(df_psO_jki_seq5_rarefied_phylum_rel, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/df_psO_jki_seq5_rarefied_phylum_rel.csv")


#For psO_jki_seq5_rarefied_annotation
psO_jki_seq5_rarefied_annotation_rel<-transform_sample_counts(psO_jki_seq5_rarefied_annotation, function(x) (x*100)/sum(x))
psO_jki_seq5_rarefied_annotation_rel
head(otu_table(psO_jki_seq5_rarefied_annotation_rel))

#### Create tables 
df_psO_jki_seq5_rarefied_annotation_rel <- data.frame(tax_table(psO_jki_seq5_rarefied_annotation_rel),otu_table(psO_jki_seq5_rarefied_annotation_rel))
write.csv(df_psO_jki_seq5_rarefied_annotation_rel, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/df_psO_jki_seq5_rarefied_annotation_rel.csv")


########################### SUBSETTING
### T1 and T2

###Subsetting samples from cleaned dataset to keep only samples that represent every case on variables  
##Variables / colors

##Using subset_samples()
psO_jki_seq5_rarefied_T1 <- subset_samples(psO_jki_seq5_rarefied, Stage=="T1")
psO_jki_seq5_rarefied_T1

psO_jki_seq5_rarefied_T2 <- subset_samples(psO_jki_seq5_rarefied, Stage=="T2")
psO_jki_seq5_rarefied_T2

###Filtering table removing ASV that do not have counts in any sample on the group subset
##Using prune_taxa()
psO_jki_seq5_rarefied_T1_filt<-prune_taxa(taxa_sums(psO_jki_seq5_rarefied_T1) > 0, psO_jki_seq5_rarefied_T1)
psO_jki_seq5_rarefied_T1_filt

psO_jki_seq5_rarefied_T2_filt<-prune_taxa(taxa_sums(psO_jki_seq5_rarefied_T2) > 0, psO_jki_seq5_rarefied_T2)
psO_jki_seq5_rarefied_T2_filt


###Agglomerating taxa per sample at the appropriated taxonomic rank
##Using tax_glom()
colnames(tax_table(psO_jki_seq5_rarefied_T1_filt))
psO_jki_seq5_rarefied_T1_filt_annotation <- tax_glom(psO_jki_seq5_rarefied_T1_filt, taxrank = "Annotation")
ntaxa(psO_jki_seq5_rarefied_T1_filt); ntaxa(psO_jki_seq5_rarefied_T1_filt_annotation)
psO_jki_seq5_rarefied_T1_filt_annotation

colnames(tax_table(psO_jki_seq5_rarefied_T2_filt))
psO_jki_seq5_rarefied_T2_filt_annotation <- tax_glom(psO_jki_seq5_rarefied_T2_filt, taxrank = "Annotation")
ntaxa(psO_jki_seq5_rarefied_T2_filt); ntaxa(psO_jki_seq5_rarefied_T2_filt_annotation)
psO_jki_seq5_rarefied_T2_filt_annotation

### Transform to relative abundance
psO_jki_seq5_rarefied_T1_filt_annotation_rel<-transform_sample_counts(psO_jki_seq5_rarefied_T1_filt_annotation, function(x) (x*100)/sum(x))
psO_jki_seq5_rarefied_T1_filt_annotation_rel
head(otu_table(psO_jki_seq5_rarefied_T1_filt_annotation_rel))

psO_jki_seq5_rarefied_T2_filt_annotation_rel<-transform_sample_counts(psO_jki_seq5_rarefied_T2_filt_annotation, function(x) (x*100)/sum(x))
psO_jki_seq5_rarefied_T2_filt_annotation_rel
head(otu_table(psO_jki_seq5_rarefied_T2_filt_annotation_rel))

### Create tables
df_psO_jki_seq5_rarefied_T1_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq5_rarefied_T1_filt_annotation_rel),otu_table(psO_jki_seq5_rarefied_T1_filt_annotation_rel))
write.csv(df_psO_jki_seq5_rarefied_T1_filt_annotation_rel, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/df_psO_jki_seq5_rarefied_T1_filt_annotation_rel.csv")

df_psO_jki_seq5_rarefied_T2_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq5_rarefied_T2_filt_annotation_rel),otu_table(psO_jki_seq5_rarefied_T2_filt_annotation_rel))
write.csv(df_psO_jki_seq5_rarefied_T2_filt_annotation_rel, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/df_psO_jki_seq5_rarefied_T2_filt_annotation_rel.csv")



########################### SUBSETTING
### D30 
### T1 and T2

###Subsetting samples from cleaned dataset to keep only samples that represent every case on variables  
##Variables / colors

##Using subset_samples()
psO_jki_seq5_rarefied_T1_D30 <- subset_samples(psO_jki_seq5_rarefied, Stage=="T1" & Depth=="D30")
psO_jki_seq5_rarefied_T1_D30

psO_jki_seq5_rarefied_T2_D30 <- subset_samples(psO_jki_seq5_rarefied, Stage=="T2" & Depth=="D30")
psO_jki_seq5_rarefied_T2_D30

###Filtering table removing ASV that do not have counts in any sample on the group subset
##Using prune_taxa()
psO_jki_seq5_rarefied_T1_D30_filt<-prune_taxa(taxa_sums(psO_jki_seq5_rarefied_T1_D30) > 0, psO_jki_seq5_rarefied_T1_D30)
psO_jki_seq5_rarefied_T1_D30_filt

psO_jki_seq5_rarefied_T2_D30_filt<-prune_taxa(taxa_sums(psO_jki_seq5_rarefied_T2_D30) > 0, psO_jki_seq5_rarefied_T2_D30)
psO_jki_seq5_rarefied_T2_D30_filt

### Create tables
df_psO_jki_seq5_rarefied_T1_D30_filt <- data.frame(tax_table(psO_jki_seq5_rarefied_T1_D30_filt),otu_table(psO_jki_seq5_rarefied_T1_D30_filt))
write.csv(df_psO_jki_seq5_rarefied_T1_D30_filt, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/df_psO_jki_seq5_rarefied_T1_D30_filt.csv")

df_psO_jki_seq5_rarefied_T2_D30_filt <- data.frame(tax_table(psO_jki_seq5_rarefied_T2_D30_filt),otu_table(psO_jki_seq5_rarefied_T2_D30_filt))
write.csv(df_psO_jki_seq5_rarefied_T2_D30_filt, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/df_psO_jki_seq5_rarefied_T2_D30_filt.csv")


###Agglomerating taxa per sample at the appropriated taxonomic rank
##Using tax_glom()
colnames(tax_table(psO_jki_seq5_rarefied_T1_D30_filt))
psO_jki_seq5_rarefied_T1_D30_filt_annotation <- tax_glom(psO_jki_seq5_rarefied_T1_D30_filt, taxrank = "Annotation")
ntaxa(psO_jki_seq5_rarefied_T1_D30_filt); ntaxa(psO_jki_seq5_rarefied_T1_D30_filt_annotation)
psO_jki_seq5_rarefied_T1_D30_filt_annotation

colnames(tax_table(psO_jki_seq5_rarefied_T2_D30_filt))
psO_jki_seq5_rarefied_T2_D30_filt_annotation <- tax_glom(psO_jki_seq5_rarefied_T2_D30_filt, taxrank = "Annotation")
ntaxa(psO_jki_seq5_rarefied_T2_D30_filt); ntaxa(psO_jki_seq5_rarefied_T2_D30_filt_annotation)
psO_jki_seq5_rarefied_T2_D30_filt_annotation

### Transform to relative abundance
psO_jki_seq5_rarefied_T1_D30_filt_annotation_rel<-transform_sample_counts(psO_jki_seq5_rarefied_T1_D30_filt_annotation, function(x) (x*100)/sum(x))
psO_jki_seq5_rarefied_T1_D30_filt_annotation_rel
head(otu_table(psO_jki_seq5_rarefied_T1_D30_filt_annotation_rel))

psO_jki_seq5_rarefied_T2_D30_filt_annotation_rel<-transform_sample_counts(psO_jki_seq5_rarefied_T2_D30_filt_annotation, function(x) (x*100)/sum(x))
psO_jki_seq5_rarefied_T2_D30_filt_annotation_rel
head(otu_table(psO_jki_seq5_rarefied_T2_D30_filt_annotation_rel))

### Create tables
df_psO_jki_seq5_rarefied_T1_D30_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq5_rarefied_T1_D30_filt_annotation_rel),otu_table(psO_jki_seq5_rarefied_T1_D30_filt_annotation_rel))
write.csv(df_psO_jki_seq5_rarefied_T1_D30_filt_annotation_rel, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/df_psO_jki_seq5_rarefied_T1_D30_filt_annotation_rel.csv")

df_psO_jki_seq5_rarefied_T2_D30_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq5_rarefied_T2_D30_filt_annotation_rel),otu_table(psO_jki_seq5_rarefied_T2_D30_filt_annotation_rel))
write.csv(df_psO_jki_seq5_rarefied_T2_D30_filt_annotation_rel, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/df_psO_jki_seq5_rarefied_T2_D30_filt_annotation_rel.csv")


               

########################### SUBSETTING
### D60 
### T1 and T2

###Subsetting samples from cleaned dataset to keep only samples that represent every case on variables  
##Variables / colors

##Using subset_samples()
psO_jki_seq5_rarefied_T1_D60 <- subset_samples(psO_jki_seq5_rarefied, Stage=="T1" & Depth=="D60")
psO_jki_seq5_rarefied_T1_D60

psO_jki_seq5_rarefied_T2_D60 <- subset_samples(psO_jki_seq5_rarefied, Stage=="T2" & Depth=="D60")
psO_jki_seq5_rarefied_T2_D60

###Filtering table removing ASV that do not have counts in any sample on the group subset
##Using prune_taxa()
psO_jki_seq5_rarefied_T1_D60_filt<-prune_taxa(taxa_sums(psO_jki_seq5_rarefied_T1_D60) > 0, psO_jki_seq5_rarefied_T1_D60)
psO_jki_seq5_rarefied_T1_D60_filt

psO_jki_seq5_rarefied_T2_D60_filt<-prune_taxa(taxa_sums(psO_jki_seq5_rarefied_T2_D60) > 0, psO_jki_seq5_rarefied_T2_D60)
psO_jki_seq5_rarefied_T2_D60_filt

### Create tables
df_psO_jki_seq5_rarefied_T1_D60_filt <- data.frame(tax_table(psO_jki_seq5_rarefied_T1_D60_filt),otu_table(psO_jki_seq5_rarefied_T1_D60_filt))
write.csv(df_psO_jki_seq5_rarefied_T1_D60_filt, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/df_psO_jki_seq5_rarefied_T1_D60_filt.csv")

df_psO_jki_seq5_rarefied_T2_D60_filt <- data.frame(tax_table(psO_jki_seq5_rarefied_T2_D60_filt),otu_table(psO_jki_seq5_rarefied_T2_D60_filt))
write.csv(df_psO_jki_seq5_rarefied_T2_D60_filt, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/df_psO_jki_seq5_rarefied_T2_D60_filt.csv")


###Agglomerating taxa per sample at the appropriated taxonomic rank
##Using tax_glom()
colnames(tax_table(psO_jki_seq5_rarefied_T1_D60_filt))
psO_jki_seq5_rarefied_T1_D60_filt_annotation <- tax_glom(psO_jki_seq5_rarefied_T1_D60_filt, taxrank = "Annotation")
ntaxa(psO_jki_seq5_rarefied_T1_D60_filt); ntaxa(psO_jki_seq5_rarefied_T1_D60_filt_annotation)
psO_jki_seq5_rarefied_T1_D60_filt_annotation

colnames(tax_table(psO_jki_seq5_rarefied_T2_D60_filt))
psO_jki_seq5_rarefied_T2_D60_filt_annotation <- tax_glom(psO_jki_seq5_rarefied_T2_D60_filt, taxrank = "Annotation")
ntaxa(psO_jki_seq5_rarefied_T2_D60_filt); ntaxa(psO_jki_seq5_rarefied_T2_D60_filt_annotation)
psO_jki_seq5_rarefied_T2_D60_filt_annotation

### Transform to relative abundance
psO_jki_seq5_rarefied_T1_D60_filt_annotation_rel<-transform_sample_counts(psO_jki_seq5_rarefied_T1_D60_filt_annotation, function(x) (x*100)/sum(x))
psO_jki_seq5_rarefied_T1_D60_filt_annotation_rel
head(otu_table(psO_jki_seq5_rarefied_T1_D60_filt_annotation_rel))

psO_jki_seq5_rarefied_T2_D60_filt_annotation_rel<-transform_sample_counts(psO_jki_seq5_rarefied_T2_D60_filt_annotation, function(x) (x*100)/sum(x))
psO_jki_seq5_rarefied_T2_D60_filt_annotation_rel
head(otu_table(psO_jki_seq5_rarefied_T2_D60_filt_annotation_rel))

### Create tables
df_psO_jki_seq5_rarefied_T1_D60_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq5_rarefied_T1_D60_filt_annotation_rel),otu_table(psO_jki_seq5_rarefied_T1_D60_filt_annotation_rel))
write.csv(df_psO_jki_seq5_rarefied_T1_D60_filt_annotation_rel, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/df_psO_jki_seq5_rarefied_T1_D60_filt_annotation_rel.csv")

df_psO_jki_seq5_rarefied_T2_D60_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq5_rarefied_T2_D60_filt_annotation_rel),otu_table(psO_jki_seq5_rarefied_T2_D60_filt_annotation_rel))
write.csv(df_psO_jki_seq5_rarefied_T2_D60_filt_annotation_rel, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/df_psO_jki_seq5_rarefied_T2_D60_filt_annotation_rel.csv")



## The end!

# psO_jki_seq2_seq4_rarefied_Go_2020 <- subset_samples(psO_jki_seq2_seq4_rarefied, Site == "Go" & Year=="Y2020")
# psO_jki_seq2_seq4_rarefied_Go_2020