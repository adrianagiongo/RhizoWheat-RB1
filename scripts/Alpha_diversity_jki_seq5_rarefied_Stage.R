#Create rarefied data from phyloseq object using vegan package or microbiome package
#Loading package
library("phyloseq")
library("vegan")
library("microbiome")
library("ggpubr")
library("ggplot2")
library("svglite")

### T1 and T2
##Calculates the alpha diversity of microbial communities
#Using dominance, richness, evenness, and diversity, all functions from microbiome package
alfa_div_psO_jki_seq5_rarefied_T1_filt_1 <- richness(psO_jki_seq5_rarefied_T1_filt, c("observed","chao1"), detection = 0)
alfa_div_psO_jki_seq5_rarefied_T1_filt_2 <- evenness(psO_jki_seq5_rarefied_T1_filt, index = 'pielou', zeroes = TRUE, detection = 0)
alfa_div_psO_jki_seq5_rarefied_T1_filt_3 <- microbiome::diversity(psO_jki_seq5_rarefied_T1_filt, index = 'shannon', zeroes = TRUE)

alfa_div_psO_jki_seq5_rarefied_T2_filt_1 <- richness(psO_jki_seq5_rarefied_T2_filt, c("observed","chao1"), detection = 0)
alfa_div_psO_jki_seq5_rarefied_T2_filt_2 <- evenness(psO_jki_seq5_rarefied_T2_filt, index = 'pielou', zeroes = TRUE, detection = 0)
alfa_div_psO_jki_seq5_rarefied_T2_filt_3 <- microbiome::diversity(psO_jki_seq5_rarefied_T2_filt, index = 'shannon', zeroes = TRUE)

##################################################
##Boxplot alpha diversity
#Prepare file using meta function
psO_jki_seq5_rarefied_T1_filt.meta <- meta(psO_jki_seq5_rarefied_T1_filt)

psO_jki_seq5_rarefied_T2_filt.meta <- meta(psO_jki_seq5_rarefied_T2_filt)

#select column from output file to be plotted and tested (only one by time)
psO_jki_seq5_rarefied_T1_filt.meta$observed <- alfa_div_psO_jki_seq5_rarefied_T1_filt_1$observed
psO_jki_seq5_rarefied_T1_filt.meta$chao1    <- alfa_div_psO_jki_seq5_rarefied_T1_filt_1$chao1
psO_jki_seq5_rarefied_T1_filt.meta$pielou   <- alfa_div_psO_jki_seq5_rarefied_T1_filt_2$pielou
psO_jki_seq5_rarefied_T1_filt.meta$shannon  <- alfa_div_psO_jki_seq5_rarefied_T1_filt_3$shannon

head(psO_jki_seq5_rarefied_T1_filt.meta)
write.csv(psO_jki_seq5_rarefied_T1_filt.meta, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/psO_jki_seq5_rarefied_T1_filt_meta.csv")

psO_jki_seq5_rarefied_T2_filt.meta$observed <- alfa_div_psO_jki_seq5_rarefied_T2_filt_1$observed
psO_jki_seq5_rarefied_T2_filt.meta$chao1    <- alfa_div_psO_jki_seq5_rarefied_T2_filt_1$chao1
psO_jki_seq5_rarefied_T2_filt.meta$pielou   <- alfa_div_psO_jki_seq5_rarefied_T2_filt_2$pielou
psO_jki_seq5_rarefied_T2_filt.meta$shannon  <- alfa_div_psO_jki_seq5_rarefied_T2_filt_3$shannon

head(psO_jki_seq5_rarefied_T2_filt.meta)
write.csv(psO_jki_seq5_rarefied_T2_filt.meta, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/psO_jki_seq5_rarefied_T2_filt_meta.csv")

#### Comparison = Rotation
#Select variable of comparison
alfa_div_comparison_rotation <- list(c("W1", "W2"))

######### T1
#Plots SHANNON
plot_psO_jki_seq5_rarefied_T1_filt_shannon_rot <- ggboxplot(psO_jki_seq5_rarefied_T1_filt.meta, 
                                                        "Depth", "shannon", fill = "Rotation", 
                                                        color = "Rotation", palette = c("#008040", "#ee4035"), 
                                                        ylab = "Shannon index", alpha = 0.8) +
  theme_bw() +  
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(6.0, 7.0), breaks = seq(6, 7, by = 0.5), label = c("6.0", "6.5", "7.0")) +
  scale_x_discrete(limits = c("D30", "D60"),  # Change these to the desired order and levels
                  labels = c("D30" = "0-30", 
                            "D60" = "30-60")) +  # Change these to desired labels
  theme(legend.position = 'right') +
  theme(legend.position = 'right') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 24, colour = "black")) +
  theme(axis.text.x = element_text(size = 24, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size = 24)) +
  theme(strip.text.x = element_text(size = 24), strip.text.y = element_text(size = 24)) +
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 24, colour = "black")) 
#stat_compare_means(method = "kruskal", label= "p", label.y = 7.0, label.x = 1, size=4) +
#stat_compare_means(method = "wilcox.test", comparisons = alfa_div_comparison_rotation, label = "p.format",  bracket.size = .2, size=3, label.y = c(6.5))
plot_psO_jki_seq5_rarefied_T1_filt_shannon_rot

ggsave("plot_psO_jki_seq5_rarefied_T1_filt_shannon_rot.png", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Alpha_div_jki_seq5/", width = 16, height = 12, units = "cm",dpi = 300)

ggsave("plot_psO_jki_seq5_rarefied_T1_filt_shannon_rot.svg", plot_psO_jki_seq5_rarefied_T1_filt_shannon_rot, device = "svg", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Alpha_div_jki_seq5/")
       
#Plots CHAO1
plot_psO_jki_seq5_rarefied_T1_filt_chao1_rot <- ggboxplot(psO_jki_seq5_rarefied_T1_filt.meta, "Depth", "chao1", fill = "Rotation", color = "Rotation", palette = c("#008040", "#ee4035"), ylab = "Chao1 index", alpha = 0.8) +
  theme_bw() + 
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(900, 2400), breaks = seq(900, 2400, by = 500), label = c("900", "1400", "1900", "2400")) +
  scale_x_discrete(limits = c("D30", "D60"),  # Change these to the desired order and levels
                   labels = c("D30" = "0-30", 
                              "D60" = "30-60")) +  # Change these to desired labels
  theme(legend.position = 'right') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) 
#stat_compare_means(method = "wilcox.test", label= "p", label.y = 2450, label.x = 1, size=4)
plot_psO_jki_seq5_rarefied_T1_filt_chao1_rot

ggsave("plot_psO_jki_seq5_rarefied_T1_filt_chao1_rot.png", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Alpha_div_jki_seq5/", width = 16, height = 12, units = "cm", dpi = 300, device = "png")

#Plots PIELOU
plot_psO_jki_seq5_rarefied_T1_filt_pielou_rot <- ggboxplot(psO_jki_seq5_rarefied_T1_filt.meta, "Depth", "pielou", fill = "Rotation", color = "Rotation", palette = c("#008040", "#ee4035"), ylab = "Pielou index", alpha = 0.8) +
  theme_bw() +  
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(0.87, 0.91), breaks = seq(0.87, 0.91, by = 0.02), label = c("0.87", "0.89", "0.91")) +
  scale_x_discrete(limits = c("D30", "D60"),  # Change these to the desired order and levels
                   labels = c("D30" = "0-30", 
                              "D60" = "30-60")) +  # Change these to desired labels
  theme(legend.position = "right") +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) 
# stat_compare_means(method = "wilcox.test", label= "p", label.y = 0.90, label.x = 1, size=4)
plot_psO_jki_seq5_rarefied_T1_filt_pielou_rot

ggsave("plot_psO_jki_seq5_rarefied_T1_filt_pielou_rot.png", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Alpha_div_jki_seq5/", width = 16, height = 12, units = "cm", dpi = 300, device = "png")


######### T2
#Plots SHANNON
plot_psO_jki_seq5_rarefied_T2_filt_shannon_rot <- ggboxplot(psO_jki_seq5_rarefied_T2_filt.meta, "Depth", "shannon", fill = "Rotation", color = "Rotation", palette = c("#008040", "#ee4035"), ylab = "Shannon index", alpha = 0.8) +
  theme_bw() +  
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(6.0, 7.0), breaks = seq(6, 7, by = 0.5), label = c("6.0", "6.5", "7.0")) +
  scale_x_discrete(limits = c("D30", "D60"),  # Change these to the desired order and levels
                   labels = c("D30" = "0-30", 
                              "D60" = "30-60")) +  # Change these to desired labels
  theme(legend.position = 'right') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 24, colour = "black")) +
  theme(axis.text.x = element_text(size = 24, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size = 24)) +
  theme(strip.text.x = element_text(size = 24), strip.text.y = element_text(size = 24)) +
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 24, colour = "black")) 
#stat_compare_means(method = "kruskal", label= "p", label.y = 7.0, label.x = 1, size=4) 
#stat_compare_means(method = "wilcox.test", comparisons = alfa_div_comparison_depth, label = "p.format",  bracket.size = .2, size=3, label.y = c(6.7, 6.8))
plot_psO_jki_seq5_rarefied_T2_filt_shannon_rot

ggsave("plot_psO_jki_seq5_rarefied_T2_filt_shannon_rot.png", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Alpha_div_jki_seq5/", width = 11, height = 12, units = "cm",dpi = 300)

#Plots CHAO1
plot_psO_jki_seq5_rarefied_T2_filt_chao1_rot <- ggboxplot(psO_jki_seq5_rarefied_T2_filt.meta, "Depth", "chao1", fill = "Rotation", color = "Rotation", palette = c("#008040", "#ee4035"), ylab = "Chao1 index", alpha = 0.8) +
  theme_bw() + 
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(900, 2400), breaks = seq(900, 2400, by = 500), label = c("900", "1400", "1900", "2400")) +
  scale_x_discrete(limits = c("D30", "D60"),  # Change these to the desired order and levels
                   labels = c("D30" = "0-30", 
                              "D60" = "30-60")) +  # Change these to desired labels
  theme(legend.position = 'right') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) 
#stat_compare_means(method = "wilcox.test", label= "p", label.y = 2450, label.x = 1, size=4)
plot_psO_jki_seq5_rarefied_T2_filt_chao1_rot

ggsave("plot_psO_jki_seq5_rarefied_T2_filt_chao1_rot.png", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Alpha_div_jki_seq5/", width = 11, height = 12, units = "cm", dpi = 300, device = "png")

#Plots PIELOU
plot_psO_jki_seq5_rarefied_T2_filt_pielou_rot <- ggboxplot(psO_jki_seq5_rarefied_T2_filt.meta, "Depth", "pielou", fill = "Rotation", color = "Rotation", palette = c("#008040", "#ee4035"), ylab = "Pielou index", alpha = 0.8) +
  theme_bw() +  
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(0.87, 0.91), breaks = seq(0.87, 0.91, by = 0.02), label = c("0.87", "0.89", "0.91")) +
  scale_x_discrete(limits = c("D30", "D60"),  # Change these to the desired order and levels
                   labels = c("D30" = "0-30", 
                              "D60" = "30-60")) +  # Change these to desired labels
  theme(legend.position = "right") +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) 
#stat_compare_means(method = "wilcox.test", label= "p", label.y = 0.90, label.x = 1, size=4)
plot_psO_jki_seq5_rarefied_T2_filt_pielou_rot

ggsave("plot_psO_jki_seq5_rarefied_T2_filt_pielou_rot.png", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Alpha_div_jki_seq5/", width = 11, height = 12, units = "cm", dpi = 300, device = "png")

## Have a good day!  : )
