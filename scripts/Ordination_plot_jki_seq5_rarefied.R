##Creating MDS plot for selected dataset
#load package
library("phyloseq")
library("vegan")
library("ggplot2")
library("dplyr")
library("microbiome")
library("ggpubr")

################################   All samples    ---------> ASV  (Figure 1B)
#set seed
set.seed(2022)

# Set color palette
Rotation_colors <- c("W1" = "#d7d3a9",  "W2" = "#74a553")

### All samples ----> MDS
MDS_bray_psO_jki_seq5_rarefied_sqr<-ordinate(psO_jki_seq5_rarefied, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq5_rarefied_T1_filt_sqr<-ordinate(psO_jki_seq5_rarefied_T1_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq5_rarefied_T2_filt_sqr<-ordinate(psO_jki_seq5_rarefied_T2_filt, "MDS","bray", autotransform=TRUE)

#Print stress data, dimensions and number of tries
head(MDS_bray_psO_jki_seq5_rarefied_sqr)
head(MDS_bray_psO_jki_seq5_rarefied_T1_filt_sqr)
head(MDS_bray_psO_jki_seq5_rarefied_T2_filt_sqr)

#Create a MDS plot 
plot_MDS_bray_psO_jki_seq5_rarefied_sqr<-plot_ordination(psO_jki_seq5_rarefied, MDS_bray_psO_jki_seq5_rarefied_sqr, type="sample",color="Rotation", shape = "Depth") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c(Rotation_colors <- c("#d7d3a9", "#74a553")))
plot_MDS_bray_psO_jki_seq5_rarefied_sqr 

ggsave("plot_MDS_bray_psO_jki_seq5_rarefied_sqr.png", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Ordination_jki_seq5/", width = 16, height = 12, units = "cm", dpi = 300, device = "png")

#Create a MDS plot 
plot_MDS_bray_psO_jki_seq5_rarefied_T1_filt_sqr<-plot_ordination(psO_jki_seq5_rarefied_T1_filt, MDS_bray_psO_jki_seq5_rarefied_T1_filt_sqr, type="sample",color="Rotation", shape = "Depth") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c(Rotation_colors <- c("#d7d3a9", "#74a553")))
plot_MDS_bray_psO_jki_seq5_rarefied_T1_filt_sqr 

ggsave("plot_MDS_bray_psO_jki_seq5_rarefied_T1_filt_sqr.png", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Ordination_jki_seq5/", width = 16, height = 12, units = "cm", dpi = 300, device = "png")


#Create a MDS plot 
plot_MDS_bray_psO_jki_seq5_rarefied_T2_filt_sqr<-plot_ordination(psO_jki_seq5_rarefied_T2_filt, MDS_bray_psO_jki_seq5_rarefied_T2_filt_sqr, type="sample",color="Rotation", shape = "Depth") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c(Depth_colors <- c("#d7d3a9", "#74a553")))
plot_MDS_bray_psO_jki_seq5_rarefied_T2_filt_sqr 

ggsave("plot_MDS_bray_psO_jki_seq5_rarefied_T2_filt_sqr.png", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Ordination_jki_seq5/", width = 16, height = 12, units = "cm", dpi = 300, device = "png")



################################   All samples    ---------> ASV    ----------> Statistics
set.seed(2022)

##Multivariate Anova for ordinate files
#Transform abundances to square root
psO_jki_seq5_rarefied_sqr <- microbiome::transform(psO_jki_seq5_rarefied, "hellinger")
psO_jki_seq5_rarefied_T1_sqr <- microbiome::transform(psO_jki_seq5_rarefied_T1_filt, "hellinger")
psO_jki_seq5_rarefied_T2_sqr <- microbiome::transform(psO_jki_seq5_rarefied_T2_filt, "hellinger")

#Convert phyloseq object to dataframe using abundances function and meta function from microbiome package
psO_jki_seq5_rarefied_sqr_abundances <- abundances(psO_jki_seq5_rarefied_sqr)
psO_jki_seq5_rarefied_meta <- meta(psO_jki_seq5_rarefied)

psO_jki_seq5_rarefied_T1_sqr_abundances <- abundances(psO_jki_seq5_rarefied_T1_sqr)
psO_jki_seq5_rarefied_T1_meta <- meta(psO_jki_seq5_rarefied_T1_filt)

psO_jki_seq5_rarefied_T2_sqr_abundances <- abundances(psO_jki_seq5_rarefied_T2_sqr)
psO_jki_seq5_rarefied_T2_meta <- meta(psO_jki_seq5_rarefied_T2_filt)

#######Permanova using adonis function from vegan and print p value
#ASVs
#by terms (it is default)
Permanova_terms_psO_jki_seq5_rarefied <- adonis2(t(psO_jki_seq5_rarefied_sqr_abundances) ~Rotation*Depth*Microhabitat*Stage, data = psO_jki_seq5_rarefied_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq5_rarefied
write.csv(Permanova_terms_psO_jki_seq5_rarefied, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/Permanova_terms_psO_jki_seq5_rarefied.csv")

Permanova_terms_psO_jki_seq5_rarefied_T1 <- adonis2(t(psO_jki_seq5_rarefied_T1_sqr_abundances) ~Microhabitat*Rotation*Depth, data = psO_jki_seq5_rarefied_T1_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq5_rarefied_T1
write.csv(Permanova_terms_psO_jki_seq5_rarefied_T1, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/Permanova_terms_psO_jki_seq5_rarefied_T1.csv")

Permanova_terms_psO_jki_seq5_rarefied_T2 <- adonis2(t(psO_jki_seq5_rarefied_T2_sqr_abundances) ~Rotation*Depth, data = psO_jki_seq5_rarefied_T2_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq5_rarefied_T2
write.csv(Permanova_terms_psO_jki_seq5_rarefied_T2, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/Permanova_terms_psO_jki_seq5_rarefied_T2.csv")


#### ANOSIM    -------> ROTATIONS
## testing of significance for the bray-curtis dissimilarity using ANOSIM
#calculate distance values between samples
dist_psO_jki_seq5_rarefied = phyloseq::distance(psO_jki_seq5_rarefied, method = "bray")

dist_psO_jki_seq5_rarefied_T1 = phyloseq::distance(psO_jki_seq5_rarefied_T1_filt, method = "bray")
dist_psO_jki_seq5_rarefied_T2 = phyloseq::distance(psO_jki_seq5_rarefied_T2_filt, method = "bray")

#create a dataframe of the metadata
metadata_psO_jki_seq5_rarefied <- data.frame(sample_data(psO_jki_seq5_rarefied))
metadata_psO_jki_seq5_rarefied_T1 <- data.frame(sample_data(psO_jki_seq5_rarefied_T1_filt))
metadata_psO_jki_seq5_rarefied_T2 <- data.frame(sample_data(psO_jki_seq5_rarefied_T2_filt))

#run Anosim using 10000 permutations
anosim_dist_psO_jki_seq5_rarefied_Rot_Stage_Mic <- anosim(dist_psO_jki_seq5_rarefied, metadata_psO_jki_seq5_rarefied$Rot_Stage_Mic, permutations = 10000)
anosim_dist_psO_jki_seq5_rarefied_T1_Rot_Stage_Mic <- anosim(dist_psO_jki_seq5_rarefied_T1, metadata_psO_jki_seq5_rarefied_T1$Rot_Stage_Mic, permutations = 10000)
anosim_dist_psO_jki_seq5_rarefied_T2_Rot_Stage_Mic <- anosim(dist_psO_jki_seq5_rarefied_T2, metadata_psO_jki_seq5_rarefied_T2$Rot_Stage_Mic, permutations = 10000)

#print results
print(anosim_dist_psO_jki_seq5_rarefied_Rot_Stage_Mic)
print(anosim_dist_psO_jki_seq5_rarefied_T1_Rot_Stage_Mic)
print(anosim_dist_psO_jki_seq5_rarefied_T2_Rot_Stage_Mic)


#### PAIRWISE ANOSIM
## testing of pairwise significance for the bray-curtis dissimilarity using ANOSIM
#Create metadata with an unique variable and run Anosim using 10000 permutations
metadata_Rot_Stage_Mic_psO_jki_seq5_rarefied <- combn(x=unique(metadata_psO_jki_seq5_rarefied$Rot_Stage_Mic), m=2)
p_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied <- c()

metadata_Rot_Stage_Mic_psO_jki_seq5_rarefied_T1 <- combn(x=unique(metadata_psO_jki_seq5_rarefied_T1$Rot_Stage_Mic), m=2)
p_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied_T1 <- c()

metadata_Rot_Stage_Mic_psO_jki_seq5_rarefied_T2 <- combn(x=unique(metadata_psO_jki_seq5_rarefied_T2$Rot_Stage_Mic), m=2)
p_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied_T2 <- c()

for (i in 1:ncol(metadata_Rot_Stage_Mic_psO_jki_seq5_rarefied)){
  ps_subs_psO_jki_seq5_rarefied <- subset_samples(psO_jki_seq5_rarefied, Rot_Stage_Mic %in% metadata_Rot_Stage_Mic_psO_jki_seq5_rarefied [,i])
  metadata_subs_psO_jki_seq5_rarefied <- data.frame(sample_data(ps_subs_psO_jki_seq5_rarefied))
  pairwise_anosim_bray_dist_psO_jki_seq5_rarefied <- anosim(phyloseq::distance(ps_subs_psO_jki_seq5_rarefied, method= "bray"), metadata_subs_psO_jki_seq5_rarefied$Rot_Stage_Mic)
  p_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied <- c(p_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied, pairwise_anosim_bray_dist_psO_jki_seq5_rarefied$signif[1])
}

for (i in 1:ncol(metadata_Rot_Stage_Mic_psO_jki_seq5_rarefied_T1)){
  ps_subs_psO_jki_seq5_rarefied_T1 <- subset_samples(psO_jki_seq5_rarefied_T1, Rot_Stage_Mic %in% metadata_Rot_Stage_Mic_psO_jki_seq5_rarefied_T1 [,i])
  metadata_subs_psO_jki_seq5_rarefied_T1 <- data.frame(sample_data(ps_subs_psO_jki_seq5_rarefied_T1))
  pairwise_anosim_bray_dist_psO_jki_seq5_rarefied_T1 <- anosim(phyloseq::distance(ps_subs_psO_jki_seq5_rarefied_T1, method= "bray"), metadata_subs_psO_jki_seq5_rarefied_T1$Rot_Stage_Mic)
  p_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied_T1 <- c(p_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied_T1, pairwise_anosim_bray_dist_psO_jki_seq5_rarefied_T1$signif[1])
}

for (i in 1:ncol(metadata_Rot_Stage_Mic_psO_jki_seq5_rarefied_T2)){
  ps_subs_psO_jki_seq5_rarefied_T2 <- subset_samples(psO_jki_seq5_rarefied_T2, Rot_Stage_Mic %in% metadata_Rot_Stage_Mic_psO_jki_seq5_rarefied_T2 [,i])
  metadata_subs_psO_jki_seq5_rarefied_T2 <- data.frame(sample_data(ps_subs_psO_jki_seq5_rarefied_T2))
  pairwise_anosim_bray_dist_psO_jki_seq5_rarefied_T2 <- anosim(phyloseq::distance(ps_subs_psO_jki_seq5_rarefied_T2, method= "bray"), metadata_subs_psO_jki_seq5_rarefied_T2$Rot_Stage_Mic)
  p_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied_T2 <- c(p_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied_T2, pairwise_anosim_bray_dist_psO_jki_seq5_rarefied_T2$signif[1])
}

#Adjust statistics values using BH
p_adj_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied <- p.adjust(p_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied, method = "BH")
p_adj_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied_T1 <- p.adjust(p_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied_T1, method = "BH")
p_adj_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied_T2 <- p.adjust(p_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied_T2, method = "BH")

#Create a table with statistics results
p_table_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied <- cbind.data.frame(t(metadata_Rot_Stage_Mic_psO_jki_seq5_rarefied), p=p_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied, p.adj=p_adj_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied)
p_table_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied_T1 <- cbind.data.frame(t(metadata_Rot_Stage_Mic_psO_jki_seq5_rarefied_T1), p=p_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied_T1, p.adj=p_adj_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied_T1)
p_table_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied_T2 <- cbind.data.frame(t(metadata_Rot_Stage_Mic_psO_jki_seq5_rarefied_T2), p=p_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied_T2, p.adj=p_adj_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied_T2)

#print results
print(pairwise_anosim_bray_dist_psO_jki_seq5_rarefied)
print(p_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied)
print(p_table_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied)

print(pairwise_anosim_bray_dist_psO_jki_seq5_rarefied_T1)
print(p_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied_T1)
print(p_table_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied_T1)

print(pairwise_anosim_bray_dist_psO_jki_seq5_rarefied_T2)
print(p_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied_T2)
print(p_table_pairwise_anosim_bray_Rot_Stage_Mic_psO_jki_seq5_rarefied_T2)

#### The end! Have fun!!
