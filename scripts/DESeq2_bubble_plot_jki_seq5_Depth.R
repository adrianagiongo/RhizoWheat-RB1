## Core - bubble plot - jki_seq5

#Load packages
library("ggplot2")
library("reshape2")

### T1_D30
####### 
#upload your data to R - exchange "Your_csv_file.csv" with the name of your csv file
bubble_sigtab_T1_D30 = read.csv("~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/input_boxplot_deseq2_T1_D30.csv", header = TRUE)
bubble_sigtab_T1_D30

#convert data frame from a "wide" format to a "long" format
bubble_sigtab_T1_D30_melt <- melt(bubble_sigtab_T1_D30, id = c("Sample", "Rotation"))
bubble_sigtab_T1_D30_melt

colours = c("#d7d3a9", "#74a553")

bubble_sigtab_T1_D30_melt$Sample <- factor(bubble_sigtab_T1_D30_melt$Sample,levels=unique(bubble_sigtab_T1_D30_melt$Sample))

plot_bubble_sigtab_T1_D30_melt = ggplot(bubble_sigtab_T1_D30_melt, aes(x = Sample, y = variable)) + 
  theme_bw() +
  geom_point(aes(size = value, fill = Rotation, stroke = 0.1), alpha = 1, shape = 21) + 
  scale_size_continuous(limits = c(0.0001, 5), range = c(0,10), breaks = c(0,0.1,0.5,1,5)) + 
  labs(x = "", y = "", size = "(%)", fill = "Rotation")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), # Remove x-axis ticks
        axis.ticks.y = element_blank(), # Remove y-axis ticks
        #axis.text.x = element_text(colour = "black", size = 12,  angle = 90, vjust = 0.3, hjust = 1), #in case we want to add sample names
        axis.text.y = element_text(colour = "black",  size = 11), 
        legend.text = element_text(size = 10, colour ="black"), 
        legend.title = element_text(size = 12), 
        panel.background = element_blank(), panel.border = element_rect(colour = "#4e4e4e", fill = NA, size = 0.5), 
        legend.position = "right", panel.grid.major.y = element_line(colour = "grey95")) +
  scale_fill_manual(values = colours, guide = guide_legend(override.aes = list(size=5)))
plot_bubble_sigtab_T1_D30_melt

ggsave("plot_bubble_sigtab_T1_D30_melt.tiff", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Bubble_plot_jki_seq5/", width = 14, height = 10, units = "cm", dpi = 300, device = "tiff")


### T2_D30
####### 
#upload your data to R - exchange "Your_csv_file.csv" with the name of your csv file
bubble_sigtab_T2_D30 = read.csv("~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/input_boxplot_deseq2_T2_D30.csv", header = TRUE)
bubble_sigtab_T2_D30

#convert data frame from a "wide" format to a "long" format
bubble_sigtab_T2_D30_melt <- melt(bubble_sigtab_T2_D30, id = c("Sample", "Rotation"))
bubble_sigtab_T2_D30_melt

colours = c("#d7d3a9", "#74a553")

bubble_sigtab_T2_D30_melt$Sample <- factor(bubble_sigtab_T2_D30_melt$Sample,levels=unique(bubble_sigtab_T2_D30_melt$Sample))

plot_bubble_sigtab_T2_D30_melt = ggplot(bubble_sigtab_T2_D30_melt, aes(x = Sample, y = variable)) + 
  theme_bw() +
  geom_point(aes(size = value, fill = Rotation, stroke = 0.1), alpha = 1, shape = 21) + 
  scale_size_continuous(limits = c(0.0001, 3), range = c(0,7), breaks = c(0,0.1,0.5,1,3)) + 
  labs(x = "", y = "", size = "(%)", fill = "Rotation")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), # Remove x-axis ticks
        axis.ticks.y = element_blank(), # Remove y-axis ticks
        axis.text.y = element_text(colour = "black",  size = 11), 
        legend.text = element_text(size = 10, colour ="black"), 
        legend.title = element_text(size = 12), 
        panel.background = element_blank(), panel.border = element_rect(colour = "#4e4e4e", fill = NA, size = 0.5), 
        legend.position = "right", panel.grid.major.y = element_line(colour = "grey95")) +
  scale_fill_manual(values = colours, guide = guide_legend(override.aes = list(size=5)))
plot_bubble_sigtab_T2_D30_melt

ggsave("plot_bubble_sigtab_T2_D30_melt.tiff", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Bubble_plot_jki_seq5/", width = 10, height = 8, units = "cm", dpi = 300, device = "tiff")



### T1_D60
####### 
#upload your data to R - exchange "Your_csv_file.csv" with the name of your csv file
bubble_sigtab_T1_D60 = read.csv("~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/input_boxplot_deseq2_T1_D60.csv", header = TRUE)
bubble_sigtab_T1_D60

#convert data frame from a "wide" format to a "long" format
bubble_sigtab_T1_D60_melt <- melt(bubble_sigtab_T1_D60, id = c("Sample", "Rotation"))
bubble_sigtab_T1_D60_melt

colours = c("#d7d3a9", "#74a553")

bubble_sigtab_T1_D60_melt$Sample <- factor(bubble_sigtab_T1_D60_melt$Sample,levels=unique(bubble_sigtab_T1_D60_melt$Sample))

plot_bubble_sigtab_T1_D60_melt = ggplot(bubble_sigtab_T1_D60_melt, aes(x = Sample, y = variable)) + 
  theme_bw() +
  geom_point(aes(size = value, fill = Rotation, stroke = 0.1), alpha = 1, shape = 21) + 
  scale_size_continuous(limits = c(0.0001, 5), range = c(0,10), breaks = c(0,0.1,0.5,1,5)) + 
  labs(x = "", y = "", size = "(%)", fill = "Rotation")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), # Remove x-axis ticks
        axis.ticks.y = element_blank(), # Remove y-axis ticks
        axis.text.y = element_text(colour = "black",  size = 11), 
        legend.text = element_text(size = 10, colour ="black"), 
        legend.title = element_text(size = 12), 
        panel.background = element_blank(), panel.border = element_rect(colour = "#4e4e4e", fill = NA, size = 0.5), 
        legend.position = "right", panel.grid.major.y = element_line(colour = "grey95")) +
  scale_fill_manual(values = colours, guide = guide_legend(override.aes = list(size=5)))
plot_bubble_sigtab_T1_D60_melt

ggsave("plot_bubble_sigtab_T1_D60_melt.tiff", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Bubble_plot_jki_seq5/", width = 14, height = 10, units = "cm", dpi = 300, device = "tiff")


### T2_D60
####### 
#upload your data to R - exchange "Your_csv_file.csv" with the name of your csv file
bubble_sigtab_T2_D60 = read.csv("~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/input_boxplot_deseq2_T2_D60.csv", header = TRUE)
bubble_sigtab_T2_D60

#convert data frame from a "wide" format to a "long" format
bubble_sigtab_T2_D60_melt <- melt(bubble_sigtab_T2_D60, id = c("Sample", "Rotation"))
bubble_sigtab_T2_D60_melt

colours = c("#d7d3a9", "#74a553")

bubble_sigtab_T2_D60_melt$Sample <- factor(bubble_sigtab_T2_D60_melt$Sample,levels=unique(bubble_sigtab_T2_D60_melt$Sample))

plot_bubble_sigtab_T2_D60_melt = ggplot(bubble_sigtab_T2_D60_melt, aes(x = Sample, y = variable)) + 
  theme_bw() +
  geom_point(aes(size = value, fill = Rotation, stroke = 0.1), alpha = 1, shape = 21) + 
  scale_size_continuous(limits = c(0.0001, 3), range = c(0,7), breaks = c(0,0.1,0.5,1,3)) + 
  labs(x = "", y = "", size = "(%)", fill = "Rotation")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), # Remove x-axis ticks
        axis.ticks.y = element_blank(), # Remove y-axis ticks
        axis.text.y = element_text(colour = "black",  size = 11), 
        legend.text = element_text(size = 10, colour ="black"), 
        legend.title = element_text(size = 12), 
        panel.background = element_blank(), panel.border = element_rect(colour = "#4e4e4e", fill = NA, size = 0.5), 
        legend.position = "right", panel.grid.major.y = element_line(colour = "grey95")) +
  scale_fill_manual(values = colours, guide = guide_legend(override.aes = list(size=5)))
plot_bubble_sigtab_T2_D60_melt

ggsave("plot_bubble_sigtab_T2_D60_melt.tiff", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Bubble_plot_jki_seq5/", width = 10, height = 7, units = "cm", dpi = 300, device = "tiff")

## Enjoy the day!

