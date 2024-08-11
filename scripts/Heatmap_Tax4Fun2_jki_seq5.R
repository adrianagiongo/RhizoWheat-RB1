#load libraries
library("ggplot2")
library("tidyr")
library("readxl")
library("dplyr")
library("RColorBrewer")


################ All taxa
### T1_5key
#load data 
N_cycle_T1_all_taxa_5key<-read_excel("~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/Heatmap_Tax4Fun2_N_cycle_all_taxa_merged_KO_T1_names_5keywords.xlsx")
N_cycle_T1_all_taxa_5key

###Preparing tidyr format table (variables in columns, observation in rows, values in cells), collecting a set of column names and place them into a single value column
N_cycle_T1_all_taxa_5key_processed<-gather(data=N_cycle_T1_all_taxa_5key,key = Tests, value = Result,aminobutyryl_CoA_ammonia_lyase:tyrosine_ammonia_lyase,-Stage,-Rot_Depth,-Sample)
head(N_cycle_T1_all_taxa_5key_processed)

###Convert data to square root 
#vt_long$Sqrt.abundance<-sqrt(vt.long$Absolute_Abundance)

#####Constructing the ggplot heatmap

### Heatmap GO_96  
#Convert variables on the x-axis to factors with desired order
#N_cycle_processed$Tests <- factor(N_cycle_processed$Tests, levels = c("Ggt", "Fusarium", "Rhizoctonia", "Cellulase", "Chitinase", "Glucanase", "Protease", "ACC", "IAA", "PO4", "Siderophore", "DAPG", "PCA", "PRND"))
N_cycle_T1_all_taxa_5key_processed$Rot_Depth <- factor(N_cycle_T1_all_taxa_5key_processed$Rot_Depth, levels = c("W1_D30", "W2_D30", "W1_D60", "W2_D60"))

#Convert variables on the x-axis to factors with desired order
#N_cycle_T1_all_taxa_5key_processed$Rot_Depth <- factor(N_cycle_T1_all_taxa_5key_processed$Rot_Depth, levels = rev(sort(unique(N_cycle_T1_all_taxa_5key_processed$Rot_Depth))))

heatmap_N_cycle_T1_all_taxa_5key_name<-ggplot(data=N_cycle_T1_all_taxa_5key_processed,mapping=aes(x=Sample, y=Tests, fill=Result))
heatmap_N_cycle_T1_all_taxa_5key_name + 
  geom_tile(colour="black",size=0.25) +
  facet_grid(~Rot_Depth, scales="free", space="free") +
  xlab(label="") +
  scale_fill_gradient(name = "Abundance",
                      low = "white", 
                      high = "#3891d6", 
                      na.value = "darkgrey")+ 
                      #limits = c(0, 0.03), # Define the scale limits here
                     # breaks = seq(0, 0.03, by = 0.01), # Define the breaks for the legend
                      #labels = c("0", "0.01", "0.02", "0.03")) + # Define the labels for the legend
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  theme(axis.ticks = element_blank()) +
  theme(axis.text.x = element_text(size = 16, colour = "black", angle = 90, hjust=1, vjust = 0.5)) +
  theme(axis.text.y = element_text(size=16, colour="black")) +
  theme(legend.title = element_text(size=16), legend.text = element_text(size = 16, colour = "black")) +
  theme(strip.text.x = element_text(size = 16, colour = "black", angle = 0, hjust=0.5)) +
  #theme(strip.text.y = element_text(size = 16, colour = "black", angle = 0)) +
    theme(strip.placement = "top", plot.title = element_text(hjust=0.5), axis.title.y = element_blank(), strip.background=element_rect(fill="#EEEEEE", color="#424242")) +
  ggtitle(label="")

ggsave("heatmap_N_cycle_T1_all_taxa_5key_name.png", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tax4Fun2_jki_seq5/", width = 45, height = 50, units = "cm",dpi = 300)

################ All taxa
### T2_5key
#load data 
N_cycle_T2_all_taxa_5key<-read_excel("~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/Heatmap_Tax4Fun2_N_cycle_all_taxa_merged_KO_T2_names_5keywords.xlsx")
N_cycle_T2_all_taxa_5key

###Preparing tidyr format table (variables in columns, observation in rows, values in cells), collecting a set of column names and place them into a single value column
N_cycle_T2_all_taxa_5key_processed<-gather(data=N_cycle_T2_all_taxa_5key,key = Tests, value = Result,aminobutyryl_CoA_ammonia_lyase:tyrosine_ammonia_lyase,-Stage,-Rot_Depth,-Sample)
head(N_cycle_T2_all_taxa_5key_processed)

###Convert data to square root 
#vt_long$Sqrt.abundance<-sqrt(vt.long$Absolute_Abundance)

#####Constructing the ggplot heatmap

### Heatmap GO_96  
#Convert variables on the x-axis to factors with desired order
#N_cycle_processed$Tests <- factor(N_cycle_processed$Tests, levels = c("Ggt", "Fusarium", "Rhizoctonia", "Cellulase", "Chitinase", "Glucanase", "Protease", "ACC", "IAA", "PO4", "Siderophore", "DAPG", "PCA", "PRND"))
N_cycle_T2_all_taxa_5key_processed$Rot_Depth <- factor(N_cycle_T2_all_taxa_5key_processed$Rot_Depth, levels = c("W1_D30", "W2_D30", "W1_D60", "W2_D60"))

#Convert variables on the x-axis to factors with desired order
#N_cycle_T2_all_taxa_5key_processed$Rot_Depth <- factor(N_cycle_T2_all_taxa_5key_processed$Rot_Depth, levels = rev(sort(unique(N_cycle_T2_all_taxa_5key_processed$Rot_Depth))))

heatmap_N_cycle_T2_all_taxa_5key_name<-ggplot(data=N_cycle_T2_all_taxa_5key_processed,mapping=aes(x=Sample, y=Tests, fill=Result))
heatmap_N_cycle_T2_all_taxa_5key_name + 
  geom_tile(colour="black",size=0.25) +
  facet_grid(~Rot_Depth, scales="free", space="free") +
  xlab(label="") +
  scale_fill_gradient(name = "Abundance",
                      low = "white", 
                      high = "#3891d6", 
                      na.value = "darkgrey")+ 
  #limits = c(0, 0.03), # Define the scale limits here
  # breaks = seq(0, 0.03, by = 0.01), # Define the breaks for the legend
  #labels = c("0", "0.01", "0.02", "0.03")) + # Define the labels for the legend
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  theme(axis.ticks = element_blank()) +
  theme(axis.text.x = element_text(size = 16, colour = "black", angle = 90, hjust=1, vjust = 0.5)) +
  theme(axis.text.y = element_text(size=16, colour="black")) +
  theme(legend.title = element_text(size=16), legend.text = element_text(size = 16, colour = "black")) +
  theme(strip.text.x = element_text(size = 16, colour = "black", angle = 0, hjust=0.5)) +
  #theme(strip.text.y = element_text(size = 16, colour = "black", angle = 0)) +
  theme(strip.placement = "top", plot.title = element_text(hjust=0.5), axis.title.y = element_blank(), strip.background=element_rect(fill="#EEEEEE", color="#424242")) +
  ggtitle(label="")

ggsave("heatmap_N_cycle_T2_all_taxa_5key_name.png", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tax4Fun2_jki_seq5/", width = 45, height = 50, units = "cm",dpi = 300)


################ Deseq2
### T1_5key
#load data 
N_cycle_T1_deseq2_5key<-read_excel("~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/Heatmap_Tax4Fun2_N_cycle_deseq2_merged_KO_T1_names_5keywords.xlsx")
N_cycle_T1_deseq2_5key

###Preparing tidyr format table (variables in columns, observation in rows, values in cells), collecting a set of column names and place them into a single value column
N_cycle_T1_deseq2_5key_processed<-gather(data=N_cycle_T1_deseq2_5key,key = Tests, value = Result,ADP_ribosyl_dinitrogen_reductase_hydrolase:two_component_system_NtrC_family_nitrogen_regulation_sensor_histidine_kinase,-Stage,-Rot_Depth,-Sample)
head(N_cycle_T1_deseq2_5key_processed)

###Convert data to square root 
#vt_long$Sqrt.abundance<-sqrt(vt.long$Absolute_Abundance)

#####Constructing the ggplot heatmap

### Heatmap GO_96  
#Convert variables on the x-axis to factors with desired order
#N_cycle_processed$Tests <- factor(N_cycle_processed$Tests, levels = c("Ggt", "Fusarium", "Rhizoctonia", "Cellulase", "Chitinase", "Glucanase", "Protease", "ACC", "IAA", "PO4", "Siderophore", "DAPG", "PCA", "PRND"))
N_cycle_T1_deseq2_5key_processed$Rot_Depth <- factor(N_cycle_T1_deseq2_5key_processed$Rot_Depth, levels = c("W1_D30", "W2_D30", "W1_D60", "W2_D60"))

#Convert variables on the x-axis to factors with desired order
#N_cycle_T1_deseq2_5key_processed$Rot_Depth <- factor(N_cycle_T1_deseq2_5key_processed$Rot_Depth, levels = rev(sort(unique(N_cycle_T1_deseq2_5key_processed$Rot_Depth))))

heatmap_N_cycle_T1_deseq2_5key_name<-ggplot(data=N_cycle_T1_deseq2_5key_processed,mapping=aes(x=Sample, y=Tests, fill=Result))
heatmap_N_cycle_T1_deseq2_5key_name + 
  geom_tile(colour="black",size=0.25) +
  facet_grid(~Rot_Depth, scales="free", space="free") +
  xlab(label="") +
  scale_fill_gradient(name = "Abundance",
                      low = "white", 
                      high = "#3891d6", 
                      na.value = "darkgrey")+ 
  #limits = c(0, 0.03), # Define the scale limits here
  # breaks = seq(0, 0.03, by = 0.01), # Define the breaks for the legend
  #labels = c("0", "0.01", "0.02", "0.03")) + # Define the labels for the legend
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  theme(axis.ticks = element_blank()) +
  theme(axis.text.x = element_text(size = 16, colour = "black", angle = 90, hjust=1, vjust = 0.5)) +
  theme(axis.text.y = element_text(size=16, colour="black")) +
  theme(legend.title = element_text(size=16), legend.text = element_text(size = 16, colour = "black")) +
  theme(strip.text.x = element_text(size = 16, colour = "black", angle = 0, hjust=0.5)) +
  #theme(strip.text.y = element_text(size = 16, colour = "black", angle = 0)) +
  theme(strip.placement = "top", plot.title = element_text(hjust=0.5), axis.title.y = element_blank(), strip.background=element_rect(fill="#EEEEEE", color="#424242")) +
  ggtitle(label="")

ggsave("heatmap_N_cycle_T1_deseq2_5key_name.png", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tax4Fun2_jki_seq5/", width = 45, height = 50, units = "cm",dpi = 300)

################ All taxa
### T2_5key
#load data 
N_cycle_T2_deseq2_5key<-read_excel("~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/Heatmap_Tax4Fun2_N_cycle_deseq2_merged_KO_T2_names_5keywords.xlsx")
N_cycle_T2_deseq2_5key

###Preparing tidyr format table (variables in columns, observation in rows, values in cells), collecting a set of column names and place them into a single value column
N_cycle_T2_deseq2_5key_processed<-gather(data=N_cycle_T2_deseq2_5key,key = Tests, value = Result,ADP_ribosyl_dinitrogen_reductase_hydrolase:two_component_system_NtrC_family_nitrogen_regulation_sensor_histidine_kinase,-Stage,-Rot_Depth,-Sample)
head(N_cycle_T2_deseq2_5key_processed)

###Convert data to square root 
#vt_long$Sqrt.abundance<-sqrt(vt.long$Absolute_Abundance)

#####Constructing the ggplot heatmap

### Heatmap GO_96  
#Convert variables on the x-axis to factors with desired order
#N_cycle_processed$Tests <- factor(N_cycle_processed$Tests, levels = c("Ggt", "Fusarium", "Rhizoctonia", "Cellulase", "Chitinase", "Glucanase", "Protease", "ACC", "IAA", "PO4", "Siderophore", "DAPG", "PCA", "PRND"))
N_cycle_T2_deseq2_5key_processed$Rot_Depth <- factor(N_cycle_T2_deseq2_5key_processed$Rot_Depth, levels = c("W1_D30", "W2_D30", "W1_D60", "W2_D60"))

#Convert variables on the x-axis to factors with desired order
#N_cycle_T2_deseq2_5key_processed$Rot_Depth <- factor(N_cycle_T2_deseq2_5key_processed$Rot_Depth, levels = rev(sort(unique(N_cycle_T2_deseq2_5key_processed$Rot_Depth))))

heatmap_N_cycle_T2_deseq2_5key_name<-ggplot(data=N_cycle_T2_deseq2_5key_processed,mapping=aes(x=Sample, y=Tests, fill=Result))
heatmap_N_cycle_T2_deseq2_5key_name + 
  geom_tile(colour="black",size=0.25) +
  facet_grid(~Rot_Depth, scales="free", space="free") +
  xlab(label="") +
  scale_fill_gradient(name = "Abundance",
                      low = "white", 
                      high = "#3891d6", 
                      na.value = "darkgrey")+ 
  #limits = c(0, 0.03), # Define the scale limits here
  # breaks = seq(0, 0.03, by = 0.01), # Define the breaks for the legend
  #labels = c("0", "0.01", "0.02", "0.03")) + # Define the labels for the legend
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  theme(axis.ticks = element_blank()) +
  theme(axis.text.x = element_text(size = 16, colour = "black", angle = 90, hjust=1, vjust = 0.5)) +
  theme(axis.text.y = element_text(size=16, colour="black")) +
  theme(legend.title = element_text(size=16), legend.text = element_text(size = 16, colour = "black")) +
  theme(strip.text.x = element_text(size = 16, colour = "black", angle = 0, hjust=0.5)) +
  #theme(strip.text.y = element_text(size = 16, colour = "black", angle = 0)) +
  theme(strip.placement = "top", plot.title = element_text(hjust=0.5), axis.title.y = element_blank(), strip.background=element_rect(fill="#EEEEEE", color="#424242")) +
  ggtitle(label="")

ggsave("heatmap_N_cycle_T2_deseq2_5key_name.png", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tax4Fun2_jki_seq5/", width = 45, height = 50, units = "cm",dpi = 300)

