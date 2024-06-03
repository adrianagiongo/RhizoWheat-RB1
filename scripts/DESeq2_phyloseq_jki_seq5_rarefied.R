##Define taxa differential abundant using DESeq2

#Loading package (to install microbiomeSeq use library("devtools")
library("phyloseq")
library("ggplot2")
library("DESeq2")
library("RColorBrewer")

### T1
##Differential analysis
#Summarize the first few entries (10) of the factor
head(sample_data(psO_jki_seq5_rarefied_T1_filt_annotation)$Rotation,10)

#Converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated
diagdds_T1 = phyloseq_to_deseq2(psO_jki_seq5_rarefied_T1_filt_annotation, ~ Rotation)

#Fitting mode and testing with Benjamini-Hochberg correction (OPTIONS - test= Wald, or t-test; and fitType = Parametric, or mean, or local)

diagdds_T1 = DESeq(diagdds_T1, test="Wald", fitType="parametric")
#diagdds_T1$Rotation <- relevel(diagdds_T1$Rotation, ref = "W1")

##Creates a table of the results of the tests stored on "diagdds", and use contrast for factors with more than 2 variables (use the reference on last position like: contrast=c("Dilution", "D3", "D0"))
res_T1 = results(diagdds_T1, contrast=c("Rotation","W2","W1"))
summary(res_T1)

##############################
#for Phyloseq plot and tables
##############################
#Indicate alpha and adjust values

#rotW1W2
alpha_DESeq2 = 0.05

sigtab_T1 = res_T1[which(res_T1$padj < alpha_DESeq2), ]

sigtab_T1 = cbind(as(sigtab_T1, "data.frame"), as(tax_table(psO_jki_seq5_rarefied_T1_filt_annotation)[rownames(sigtab_T1), ], "matrix"), as(otu_table(psO_jki_seq5_rarefied_T1_filt_annotation_rel)[rownames(sigtab_T1), ], "matrix"))
head(sigtab_T1)
dim(sigtab_T1)

write.csv(sigtab_T1, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/sigtab_T1.csv")


### T2
##Differential analysis
#Summarize the first few entries (10) of the factor
head(sample_data(psO_jki_seq5_rarefied_T2_filt_annotation)$Rotation,10)

#Converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated
diagdds_T2 = phyloseq_to_deseq2(psO_jki_seq5_rarefied_T2_filt_annotation, ~ Rotation)

#Fitting mode and testing with Benjamini-Hochberg correction (OPTIONS - test= Wald, or t-test; and fitType = Parametric, or mean, or local)

diagdds_T2 = DESeq(diagdds_T2, test="Wald", fitType="parametric")
#diagdds_T2$Rotation <- relevel(diagdds_T2$Rotation, ref = "W1")

##Creates a table of the results of the tests stored on "diagdds", and use contrast for factors with more than 2 variables (use the reference on last position like: contrast=c("Dilution", "D3", "D0"))
res_T2 = results(diagdds_T2, contrast=c("Rotation","W2","W1"))
summary(res_T2)

##############################
#for Phyloseq plot and tables
##############################
#Indicate alpha and adjust values

#rotW1W2
alpha_DESeq2 = 0.05

sigtab_T2 = res_T2[which(res_T2$padj < alpha_DESeq2), ]

sigtab_T2 = cbind(as(sigtab_T2, "data.frame"), as(tax_table(psO_jki_seq5_rarefied_T2_filt_annotation)[rownames(sigtab_T2), ], "matrix"), as(otu_table(psO_jki_seq5_rarefied_T2_filt_annotation_rel)[rownames(sigtab_T2), ], "matrix"))
head(sigtab_T2)
dim(sigtab_T2)

write.csv(sigtab_T2, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/sigtab_T2.csv")


###################
# Subsetting Stage and Depth
### T1_D30
##Differential analysis
#Summarize the first few entries (10) of the factor
head(sample_data(psO_jki_seq5_rarefied_T1_D30_filt_annotation)$Rotation,10)

#Converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated
diagdds_T1_D30 = phyloseq_to_deseq2(psO_jki_seq5_rarefied_T1_D30_filt_annotation, ~ Rotation)

#Fitting mode and testing with Benjamini-Hochberg correction (OPTIONS - test= Wald, or t-test; and fitType = Parametric, or mean, or local)

diagdds_T1_D30 = DESeq(diagdds_T1_D30, test="Wald", fitType="parametric")
#diagdds_T1_D30$Rotation <- relevel(diagdds_T1_D30$Rotation, ref = "W1")

##Creates a table of the results of the tests stored on "diagdds", and use contrast for factors with more than 2 variables (use the reference on last position like: contrast=c("Dilution", "D3", "D0"))
res_T1_D30 = results(diagdds_T1_D30, contrast=c("Rotation","W2","W1"))
summary(res_T1_D30)

##############################
#for Phyloseq plot and tables
##############################
#Indicate alpha and adjust values

#rotW1W2
alpha_DESeq2 = 0.05

sigtab_T1_D30 = res_T1_D30[which(res_T1_D30$padj < alpha_DESeq2), ]

sigtab_T1_D30 = cbind(as(sigtab_T1_D30, "data.frame"), as(tax_table(psO_jki_seq5_rarefied_T1_D30_filt_annotation)[rownames(sigtab_T1_D30), ], "matrix"), as(otu_table(psO_jki_seq5_rarefied_T1_D30_filt_annotation_rel)[rownames(sigtab_T1_D30), ], "matrix"))
head(sigtab_T1_D30)
dim(sigtab_T1_D30)

write.csv(sigtab_T1_D30, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/sigtab_T1_D30.csv")


### T2_D30
##Differential analysis
#Summarize the first few entries (10) of the factor
head(sample_data(psO_jki_seq5_rarefied_T2_D30_filt_annotation)$Rotation,10)

#Converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated
diagdds_T2_D30 = phyloseq_to_deseq2(psO_jki_seq5_rarefied_T2_D30_filt_annotation, ~ Rotation)

#Fitting mode and testing with Benjamini-Hochberg correction (OPTIONS - test= Wald, or t-test; and fitType = Parametric, or mean, or local)

diagdds_T2_D30 = DESeq(diagdds_T2_D30, test="Wald", fitType="parametric")
#diagdds_T2_D30$Rotation <- relevel(diagdds_T2_D30$Rotation, ref = "W1")

##Creates a table of the results of the tests stored on "diagdds", and use contrast for factors with more than 2 variables (use the reference on last position like: contrast=c("Dilution", "D3", "D0"))
res_T2_D30 = results(diagdds_T2_D30, contrast=c("Rotation","W2","W1"))
summary(res_T2_D30)

##############################
#for Phyloseq plot and tables
##############################
#Indicate alpha and adjust values

#rotW1W2
alpha_DESeq2 = 0.05

sigtab_T2_D30 = res_T2_D30[which(res_T2_D30$padj < alpha_DESeq2), ]

sigtab_T2_D30 = cbind(as(sigtab_T2_D30, "data.frame"), as(tax_table(psO_jki_seq5_rarefied_T2_D30_filt_annotation)[rownames(sigtab_T2_D30), ], "matrix"), as(otu_table(psO_jki_seq5_rarefied_T2_D30_filt_annotation_rel)[rownames(sigtab_T2_D30), ], "matrix"))
head(sigtab_T2_D30)
dim(sigtab_T2_D30)

write.csv(sigtab_T2_D30, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/sigtab_T2_D30.csv")



###################
# Subsetting Stage and Depth
### T1_D60
##Differential analysis
#Summarize the first few entries (10) of the factor
head(sample_data(psO_jki_seq5_rarefied_T1_D60_filt_annotation)$Rotation,10)

#Converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated
diagdds_T1_D60 = phyloseq_to_deseq2(psO_jki_seq5_rarefied_T1_D60_filt_annotation, ~ Rotation)

#Fitting mode and testing with Benjamini-Hochberg correction (OPTIONS - test= Wald, or t-test; and fitType = Parametric, or mean, or local)

diagdds_T1_D60 = DESeq(diagdds_T1_D60, test="Wald", fitType="parametric")
#diagdds_T1_D60$Rotation <- relevel(diagdds_T1_D60$Rotation, ref = "W1")

##Creates a table of the results of the tests stored on "diagdds", and use contrast for factors with more than 2 variables (use the reference on last position like: contrast=c("Dilution", "D3", "D0"))
res_T1_D60 = results(diagdds_T1_D60, contrast=c("Rotation","W2","W1"))
summary(res_T1_D60)

##############################
#for Phyloseq plot and tables
##############################
#Indicate alpha and adjust values

#rotW1W2
alpha_DESeq2 = 0.05

sigtab_T1_D60 = res_T1_D60[which(res_T1_D60$padj < alpha_DESeq2), ]

sigtab_T1_D60 = cbind(as(sigtab_T1_D60, "data.frame"), as(tax_table(psO_jki_seq5_rarefied_T1_D60_filt_annotation)[rownames(sigtab_T1_D60), ], "matrix"), as(otu_table(psO_jki_seq5_rarefied_T1_D60_filt_annotation_rel)[rownames(sigtab_T1_D60), ], "matrix"))
head(sigtab_T1_D60)
dim(sigtab_T1_D60)

write.csv(sigtab_T1_D60, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/sigtab_T1_D60.csv")


### T2_D60
##Differential analysis
#Summarize the first few entries (10) of the factor
head(sample_data(psO_jki_seq5_rarefied_T2_D60_filt_annotation)$Rotation,10)

#Converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated
diagdds_T2_D60 = phyloseq_to_deseq2(psO_jki_seq5_rarefied_T2_D60_filt_annotation, ~ Rotation)

#Fitting mode and testing with Benjamini-Hochberg correction (OPTIONS - test= Wald, or t-test; and fitType = Parametric, or mean, or local)

diagdds_T2_D60 = DESeq(diagdds_T2_D60, test="Wald", fitType="parametric")
#diagdds_T2_D60$Rotation <- relevel(diagdds_T2_D60$Rotation, ref = "W1")

##Creates a table of the results of the tests stored on "diagdds", and use contrast for factors with more than 2 variables (use the reference on last position like: contrast=c("Dilution", "D3", "D0"))
res_T2_D60 = results(diagdds_T2_D60, contrast=c("Rotation","W2","W1"))
summary(res_T2_D60)

##############################
#for Phyloseq plot and tables
##############################
#Indicate alpha and adjust values

#rotW1W2
alpha_DESeq2 = 0.05

sigtab_T2_D60 = res_T2_D60[which(res_T2_D60$padj < alpha_DESeq2), ]

sigtab_T2_D60 = cbind(as(sigtab_T2_D60, "data.frame"), as(tax_table(psO_jki_seq5_rarefied_T2_D60_filt_annotation)[rownames(sigtab_T2_D60), ], "matrix"), as(otu_table(psO_jki_seq5_rarefied_T2_D60_filt_annotation_rel)[rownames(sigtab_T2_D60), ], "matrix"))
head(sigtab_T2_D60)
dim(sigtab_T2_D60)

write.csv(sigtab_T2_D60, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/sigtab_T2_D60.csv")

##The end!

