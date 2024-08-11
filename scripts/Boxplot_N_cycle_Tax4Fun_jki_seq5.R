#Load packages
library("readxl")
library("ggplot2")
library("ggpubr")

########### Nitrogen ---> all dataset included
All4groups_all_taxa_5key_results <- read.csv ("~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/nitrogen_all_taxa_5key.csv")

##check for the distribution 
#"shapiro.test" function to test the Null hypothesis, and "qqnorm" function to qq plot
# If alpha result is lower than the alpha value chosen (p<0.05),
# then the null hypotesis (the population is normally distributed) is rejected
hist(All4groups_all_taxa_5key_results$Nitrogen)
shapiro.test(All4groups_all_taxa_5key_results$Nitrogen)
qqnorm(All4groups_all_taxa_5key_results$Nitrogen)

################ 

#Select variable of comparison
comparison_time <- list(c("W1_D30", "W2_D30"), c("W1_D60", "W2_D60"))
comparison_depth <- list(c("W1_D30", "W1_D60"), c("W2_D30", "W2_D60"))

boxplot_All4groups_all_taxa_5key <- ggboxplot(All4groups_all_taxa_5key_results, "Treatment", "Nitrogen", fill = "#008040",
                                   color = "Treatment", palette = c("#008040","#008040","#008040","#008040"),         
                                   ylab = "Nitrogen",
                                   bxp.errorbar = TRUE, bxp.errorbar.width = 0.2) +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Stage, scales = "free", space="free_x") +
  scale_y_continuous(limits = c(0.7, 1.3), breaks = seq(0.7, 1.3, by = 0.1), label = c("0.7", "0.8", "0.9", "1.0", "1.1", "1.2", "1.3")) +
  stat_compare_means(test = "kruskal.test", comparisons = comparison_time, label= "p", bracket.size = .3, size=4, label.y = c(1.1, 1.15)) +
  stat_compare_means(test = "kruskal.test", comparisons = comparison_depth, label= "p", bracket.size = .3, size=4, label.y = c(1.2, 1.2)) +
  #stat_compare_means(aes(label = after_stat(p.format)), method = "wilcox.test", ref.group = "T1_D30", label.y = 0.031) +
  #stat_compare_means(label.y = 1) + 
  ylab("Relative abundance (%)") +
  theme(legend.position = 'none') +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black"))
boxplot_All4groups_all_taxa_5key

ggsave("boxplot_All4groups_all_taxa_5key.png", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tax4Fun2_jki_seq5/", width = 20, height = 16, units = "cm",dpi = 200)



########### Nitrogen ---> all dataset included
All4groups_deseq2_5key_results <- read.csv ("~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/nitrogen_deseq2_5key.csv")

##check for the distribution 
#"shapiro.test" function to test the Null hypothesis, and "qqnorm" function to qq plot
# If alpha result is lower than the alpha value chosen (p<0.05),
# then the null hypotesis (the population is normally distributed) is rejected
hist(All4groups_deseq2_5key_results$Nitrogen)
shapiro.test(All4groups_deseq2_5key_results$Nitrogen)
qqnorm(All4groups_deseq2_5key_results$Nitrogen)

################ 

#Select variable of comparison
comparison_time <- list(c("W1_D30", "W2_D30"), c("W1_D60", "W2_D60"))
comparison_depth <- list(c("W1_D30", "W1_D60"), c("W2_D30", "W2_D60"))

boxplot_All4groups_deseq2_5key <- ggboxplot(All4groups_deseq2_5key_results, "Treatment", "Nitrogen", fill = "#3891d6",
                                            color = "Treatment", palette = c("#3891d6", "#3891d6","#3891d6","#3891d6"),         
                                            ylab = "Nitrogen",
                                            bxp.errorbar = TRUE, bxp.errorbar.width = 0.2) +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Stage, scales = "free", space="free_x") +
  scale_y_continuous(limits = c(0.7, 1.3), breaks = seq(0.7, 1.3, by = 0.1), label = c("0.7", "0.8", "0.9", "1.0", "1.1", "1.2", "1.3")) +
  stat_compare_means(test = "kruskal.test", comparisons = comparison_time, label= "p", bracket.size = .3, size=4, label.y = c(1.1, 1.15)) +
  stat_compare_means(test = "kruskal.test", comparisons = comparison_depth, label= "p", bracket.size = .3, size=4, label.y = c(1.2, 1.2)) +
  #stat_compare_means(aes(label = after_stat(p.format)), method = "wilcox.test", ref.group = "T1_D30", label.y = 0.031) +
  #stat_compare_means(label.y = 1) + 
  ylab("Relative abundance (%)") +
  theme(legend.position = 'none') +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black"))
boxplot_All4groups_deseq2_5key

ggsave("boxplot_All4groups_deseq2_5key.png", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tax4Fun2_jki_seq5/", width = 20, height = 16, units = "cm",dpi = 200)





