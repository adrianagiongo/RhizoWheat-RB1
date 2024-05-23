#https://jkzorz.github.io/2019/07/02/Indicator-species-analysis.html
#https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/isa/
  
library("phyloseq")
library("vegan")
library("microbiome")
library("ggpubr")
library("ggplot2")
library("tidyr")
library("indicspecies")

#### 2020
##Ki  --> BS
#######################################
psO_jki_seq5_rarefied_T1_filt_annotation_rel

df_psO_jki_seq5_rarefied_T1_filt_annotation_rel
str(df_psO_jki_seq5_rarefied_T1_filt_annotation_rel)

###Only annotation name
df_psO_jki_seq5_rarefied_T1_filt_annotation_rel_only <- df_psO_jki_seq5_rarefied_T1_filt_annotation_rel[,8:ncol(df_psO_jki_seq5_rarefied_T1_filt_annotation_rel)]
str(df_psO_jki_seq5_rarefied_T1_filt_annotation_rel_only)

df_psO_jki_seq5_rarefied_T1_filt_annotation_rel_only_t <- t(df_psO_jki_seq5_rarefied_T1_filt_annotation_rel_only)
str(df_psO_jki_seq5_rarefied_T1_filt_annotation_rel_only_t)

### Create table only annotation 
write.csv(df_psO_jki_seq5_rarefied_T1_filt_annotation_rel_only, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/df_psO_jki_seq5_rarefied_T1_filt_annotation_rel_only.csv")
write.csv(df_psO_jki_seq5_rarefied_T1_filt_annotation_rel_only_t, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/df_psO_jki_seq5_rarefied_T1_filt_annotation_rel_only_t.csv")

### METADATA file   ----> Rotation
#Prepare file with meta data using meta function
psO_jki_seq5_rarefied_T1_filt_annotation_rel_met.meta <- meta(psO_jki_seq5_rarefied_T1_filt_annotation_rel)
head(psO_jki_seq5_rarefied_T1_filt_annotation_rel_met.meta)
psO_jki_seq5_rarefied_T1_filt_annotation_rel_met.meta_Rotation <- subset(psO_jki_seq5_rarefied_T1_filt_annotation_rel_met.meta, 
                                                                                select = -c(Sample_name, Location, Stage, Microhabitat, Depth, Replicate,
                                                                                            Rot_Stage, Rot_Mic, Rot_Stage_Mic))
head(psO_jki_seq5_rarefied_T1_filt_annotation_rel_met.meta_Rotation)

#Combine matrix in one
combined_jki_seq5_rarefied_T1_filt_annotation_rel<-cbind(psO_jki_seq5_rarefied_T1_filt_annotation_rel_met.meta_Rotation, df_psO_jki_seq5_rarefied_T1_filt_annotation_rel_only_t)
combined_jki_seq5_rarefied_T1_filt_annotation_rel[1:8,1:5]

# Check rows and columns and structure
colnames(combined_jki_seq5_rarefied_T1_filt_annotation_rel)
rownames(combined_jki_seq5_rarefied_T1_filt_annotation_rel)
str(combined_jki_seq5_rarefied_T1_filt_annotation_rel)

#Create a data frame with all my abundance data, and a vector that contains the information from the “Site_soil” column
combined_jki_seq5_rarefied_T1_filt_annotation_rel_Rotation = combined_jki_seq5_rarefied_T1_filt_annotation_rel[,3:ncol(combined_jki_seq5_rarefied_T1_filt_annotation_rel)]
str(combined_jki_seq5_rarefied_T1_filt_annotation_rel_Rotation)
Rotation = combined_jki_seq5_rarefied_T1_filt_annotation_rel$Rotation

# Indicator species analysis (with site group combinations) 
indval_combined_jki_seq5_rarefied_T1_filt_annotation_rel_Rotation <- multipatt(combined_jki_seq5_rarefied_T1_filt_annotation_rel_Rotation, Rotation, func = "IndVal.g", control = how(nperm=999)) 
summary(indval_combined_jki_seq5_rarefied_T1_filt_annotation_rel_Rotation, indvalcomp = TRUE)

# Capture the output
indval_combined_jki_seq5_rarefied_T1_filt_annotation_rel_Rotation <- capture.output({
  summary(indval_combined_jki_seq5_rarefied_T1_filt_annotation_rel_Rotation, indvalcomp = TRUE)
})

# Write the output to a file
write.csv(indval_combined_jki_seq5_rarefied_T1_filt_annotation_rel_Rotation, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/indval_combined_jki_seq5_rarefied_T1_filt_annotation_rel_Rotation.csv")

