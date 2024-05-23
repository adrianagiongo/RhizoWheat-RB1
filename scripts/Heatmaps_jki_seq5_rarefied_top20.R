##Creating heatmap for selected dataset
rmarkdown::render("~/Documents/R_analysis/jki_seq5/scripts_jki_seq5/Heatmaps_jki_seq5.R")

##load libraries
library("phyloseq")
library("microbiome")
library("ggplot2")
library("RColorBrewer")
library("gridExtra")

## Heatmaps with annotations with TOP 10
##for psO_jki_seq2_Go_filt_annotation_rel
#List top abundant taxa
psO_jki_seq5_rarefied_T1_filt_annotation_rel_top10 <- top_taxa(psO_jki_seq5_rarefied_T1_filt_annotation_rel, n = 10)
psO_jki_seq5_rarefied_T1_filt_annotation_rel_top10

#Subset top abundant taxa
psO_jki_seq5_rarefied_T1_filt_annotation_rel_top10_taxa <- prune_taxa(psO_jki_seq5_rarefied_T1_filt_annotation_rel_top10, psO_jki_seq5_rarefied_T1_filt_annotation_rel)

#Create tables
df_psO_jki_seq5_rarefied_T1_filt_annotation_rel_top10_taxa <- data.frame(tax_table(psO_jki_seq5_rarefied_T1_filt_annotation_rel_top10_taxa),otu_table(psO_jki_seq5_rarefied_T1_filt_annotation_rel_top10_taxa))
write.csv(df_psO_jki_seq5_rarefied_T1_filt_annotation_rel_top10_taxa, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/df_psO_jki_seq5_rarefied_T1_filt_annotation_rel_top10_taxa.csv")

##Using plot_Heatmap()
###Heatmap for Go
heatmap_psO_jki_seq5_rarefied_T1_filt_annotation_rel_top10_taxa = plot_heatmap(psO_jki_seq5_rarefied_T1_filt_annotation_rel_top10_taxa, sample.order = "Rot_Stage",sample.label = "Sample_name", taxa.label = "Annotation",taxa.order="Annotation",low="white", high="#1d5b65", na.value="white", trans = NULL) +
  theme_bw() +
  #facet_grid(~Microhabitat, scales = "free", space = "free") +
  labs(fill="Relative Abundance") +
  theme(axis.text.x = element_text(size = 24, angle = 90, vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size=24, face = "italic")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.text = element_text(size=24), legend.title = element_blank()) +
  theme(strip.text.x = element_text(size=24)) +
  #scale_fill_gradientn (limits = c(0, 7), breaks = seq(0, 7, by = 3), label = c("0", "3", "6"), colours = c("white", "#1d5b65")) +
  #scale_y_discrete("Annotation", labels = c('Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium' = "ANPR group")) +
  geom_tile(colour="gray",size=0.15)
heatmap_psO_jki_seq5_rarefied_T1_filt_annotation_rel_top10_taxa

ggsave("heatmap_psO_jki_seq5_rarefied_T1_filt_annotation_rel_top10_taxa.png", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Heatmap_jki_seq5/", width = 18.5, height = 10, units = "cm",dpi = 300)

## Heatmaps with annotations with TOP 10
##for psO_jki_seq2_Go_filt_annotation_rel
#List top abundant taxa
psO_jki_seq5_rarefied_T2_filt_annotation_rel_top10 <- top_taxa(psO_jki_seq5_rarefied_T2_filt_annotation_rel, n = 10)
psO_jki_seq5_rarefied_T2_filt_annotation_rel_top10

#Subset top abundant taxa
psO_jki_seq5_rarefied_T2_filt_annotation_rel_top10_taxa <- prune_taxa(psO_jki_seq5_rarefied_T2_filt_annotation_rel_top10, psO_jki_seq5_rarefied_T2_filt_annotation_rel)

#Create tables
df_psO_jki_seq5_rarefied_T2_filt_annotation_rel_top10_taxa <- data.frame(tax_table(psO_jki_seq5_rarefied_T2_filt_annotation_rel_top10_taxa),otu_table(psO_jki_seq5_rarefied_T2_filt_annotation_rel_top10_taxa))
write.csv(df_psO_jki_seq5_rarefied_T2_filt_annotation_rel_top10_taxa, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/df_psO_jki_seq5_rarefied_T2_filt_annotation_rel_top10_taxa.csv")

##Using plot_Heatmap()
###Heatmap for Go
heatmap_psO_jki_seq5_rarefied_T2_filt_annotation_rel_top10_taxa = plot_heatmap(psO_jki_seq5_rarefied_T2_filt_annotation_rel_top10_taxa, sample.order = "Rot_Stage",sample.label = "Sample_name", taxa.label = "Annotation",taxa.order="Annotation",low="white", high="#1d5b65", na.value="white", trans = NULL) +
  theme_bw() +
  #facet_grid(~Microhabitat, scales = "free", space = "free") +
  labs(fill="Relative Abundance") +
  theme(axis.text.x = element_text(size = 24, angle = 90, vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size=24, face = "italic")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.text = element_text(size=24), legend.title = element_blank()) +
  theme(strip.text.x = element_text(size=24)) +
  #scale_fill_gradientn (limits = c(0, 7), breaks = seq(0, 7, by = 3), label = c("0", "3", "6"), colours = c("white", "#1d5b65")) +
  #scale_y_discrete("Annotation", labels = c('Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium' = "ANPR group")) +
  geom_tile(colour="gray",size=0.15)
heatmap_psO_jki_seq5_rarefied_T2_filt_annotation_rel_top10_taxa

ggsave("heatmap_psO_jki_seq5_rarefied_T2_filt_annotation_rel_top10_taxa.png", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Heatmap_jki_seq5/", width = 18.5, height = 10, units = "cm",dpi = 300)


## Enjoy the day!  : )

