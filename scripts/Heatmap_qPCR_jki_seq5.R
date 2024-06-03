#load libraries
library("ggplot2")
library("tidyr")
library("readxl")
library("dplyr")
library("RColorBrewer")

### qPCR_T1_D30
#load data 
qPCR_T1_D30<-read_excel("~/Documents/R_analysis/jki_seq5/data_jki_seq5/qPCR_bacteria_T1_D30.xlsx")
qPCR_T1_D30

###Preparing tidyr format table (variables in columns, observation in rows, values in cells), collecting a set of column names and place them into a single value column
qPCR_T1_D30_processed<-gather(data=qPCR_T1_D30,key = Tests, value = Result,amoA:nxrAB)
head(qPCR_T1_D30_processed)

###Convert data to square root 
#vt_long$Sqrt.abundance<-sqrt(vt.long$Absolute_Abundance)

#####Constructing the ggplot heatmap

### Heatmap GO_96  
#Convert variables on the x-axis to factors with desired order
qPCR_T1_D30_processed$Tests <- factor(qPCR_T1_D30_processed$Tests, levels = c("amoA", "nifH", "nirS", "nosZ", "napA", "narG", "nasA", "nirB", "nirK", "norB", "nrfA", "nxrAB"))
qPCR_T1_D30_processed$Genus_id <- factor(qPCR_T1_D30_processed$Genus_id, levels = c("Moraxella",
                                                                                    "Neisseria",
                                                                                    "Corynebacterium",
                                                                                    "Haemophilus",
                                                                                    "WWE3",
                                                                                    "Brevibacillus",
                                                                                    "Noviherbaspirillum",
                                                                                    "Porphyrobacter",
                                                                                    "Phenylobacterium",
                                                                                    "Sphingomonas",
                                                                                    "Xanthobacteraceae",
                                                                                    "Tumebacillus"))

#qPCR_T1_D30_processed$Compartment <- factor(qPCR_T1_D30_processed$Compartment, levels = c("BS", "RH", "RP"))

#Convert variables on the x-axis to factors with desired order
qPCR_T1_D30_processed$Genus_id <- factor(qPCR_T1_D30_processed$Genus_id, levels = rev(sort(unique(qPCR_T1_D30_processed$Genus_id))))

heatmap_qPCR_T1_D30<-ggplot(data=qPCR_T1_D30_processed,mapping=aes(x=Tests, y=Genus_id, fill=Result))
heatmap_qPCR_T1_D30 + 
  geom_tile(colour="darkgrey",size=0.25) +
  facet_grid(~Tests, scales="free", space="free") +
  theme(legend.position = "none") +
  xlab(label="") +
  scale_fill_gradient(name="Index",low="white", high="black", na.value="lightgrey") +
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  theme(axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size=16, colour="black", face = "italic")) +
  #theme(legend.title = element_text(size=16), legend.text = element_text(size = 16, colour = "black")) +
  theme(legend.title = element_blank()) +
  theme(strip.text.x = element_text(size = 16, colour = "black", angle = 90, hjust=0)) +
  theme(strip.text.y = element_text(size = 16, colour = "black", angle = 0)) +
  theme(strip.placement = "top", plot.title = element_text(hjust=0.5), axis.title.y = element_blank(), strip.background=element_rect(fill="#EEEEEE", color="#4e4e4e")) +
  ggtitle(label="")

ggsave("heatmap_qPCR_T1_D30.png", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Heatmap_jki_seq5/", width = 20, height = 14, units = "cm",dpi = 300)



### qPCR_T2_D30
#load data 
qPCR_T2_D30<-read_excel("~/Documents/R_analysis/jki_seq5/data_jki_seq5/qPCR_bacteria_T2_D30.xlsx")
qPCR_T2_D30

###Preparing tidyr format table (variables in columns, observation in rows, values in cells), collecting a set of column names and place them into a single value column
qPCR_T2_D30_processed<-gather(data=qPCR_T2_D30,key = Tests, value = Result,amoA:nxrAB)
head(qPCR_T2_D30_processed)

###Convert data to square root 
#vt_long$Sqrt.abundance<-sqrt(vt.long$Absolute_Abundance)

#####Constructing the ggplot heatmap

### Heatmap GO_96  
#Convert variables on the x-axis to factors with desired order
qPCR_T2_D30_processed$Tests <- factor(qPCR_T2_D30_processed$Tests, levels = c("amoA", "nifH", "nirS", "nosZ", "napA", "narG", "nasA", "nirB", "nirK", "norB", "nrfA", "nxrAB"))
qPCR_T2_D30_processed$Genus_id <- factor(qPCR_T2_D30_processed$Genus_id, levels = c("Shinella",
                                                                                    "Subgroup_10",
                                                                                    "ANPR",
                                                                                    "Acinetobacter",
                                                                                    "Devosiaceae"))

#qPCR_T2_D30_processed$Compartment <- factor(qPCR_T2_D30_processed$Compartment, levels = c("BS", "RH", "RP"))

#Convert variables on the x-axis to factors with desired order
qPCR_T2_D30_processed$Genus_id <- factor(qPCR_T2_D30_processed$Genus_id, levels = rev(sort(unique(qPCR_T2_D30_processed$Genus_id))))

heatmap_qPCR_T2_D30<-ggplot(data=qPCR_T2_D30_processed,mapping=aes(x=Tests, y=Genus_id, fill=Result))
heatmap_qPCR_T2_D30 + 
  geom_tile(colour="darkgrey",size=0.25) +
  facet_grid(~Tests, scales="free", space="free") +
  theme(legend.position = "none") +
  xlab(label="") +
  scale_fill_gradient(name="Index",low="white", high="black", na.value="lightgrey") +
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  theme(axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size=16, colour="black", face = "italic")) +
  #theme(legend.title = element_text(size=16), legend.text = element_text(size = 16, colour = "black")) +
  theme(legend.title = element_blank()) +
  theme(strip.text.x = element_text(size = 16, colour = "black", angle = 90, hjust=0)) +
  theme(strip.text.y = element_text(size = 16, colour = "black", angle = 0)) +
  theme(strip.placement = "top", plot.title = element_text(hjust=0.5), axis.title.y = element_blank(), strip.background=element_rect(fill="#EEEEEE", color="#424242")) +
  ggtitle(label="")

ggsave("heatmap_qPCR_T2_D30.png", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Heatmap_jki_seq5/", width = 20, height = 10, units = "cm",dpi = 300)



### qPCR_T1_D60
#load data 
qPCR_T1_D60<-read_excel("~/Documents/R_analysis/jki_seq5/data_jki_seq5/qPCR_bacteria_T1_D60.xlsx")
qPCR_T1_D60

###Preparing tidyr format table (variables in columns, observation in rows, values in cells), collecting a set of column names and place them into a single value column
qPCR_T1_D60_processed<-gather(data=qPCR_T1_D60,key = Tests, value = Result,amoA:nxrAB)
head(qPCR_T1_D60_processed)

###Convert data to square root 
#vt_long$Sqrt.abundance<-sqrt(vt.long$Absolute_Abundance)

#####Constructing the ggplot heatmap

### Heatmap GO_96  
#Convert variables on the x-axis to factors with desired order
qPCR_T1_D60_processed$Tests <- factor(qPCR_T1_D60_processed$Tests, levels = c("amoA", "nifH", "nirS", "nosZ", "napA", "narG", "nasA", "nirB", "nirK", "norB", "nrfA", "nxrAB"))
qPCR_T1_D60_processed$Genus_id <- factor(qPCR_T1_D60_processed$Genus_id, levels = c('Arthrospira',
                                                                                    'Blastomonas',
                                                                                    'Ktedonobacteraceae',
                                                                                    'Candidatus_Pacebacteria',
                                                                                    'Candidatus_Doudnabacteria',
                                                                                    'Candidatus_Liptonbacteria',
                                                                                    'Parcubacteria',
                                                                                    'KD3_93',
                                                                                    'Porphyrobacter',
                                                                                    'Phenylobacterium',
                                                                                    'Candidatus_Udaeobacter',
                                                                                    'Xanthobacteraceae',
                                                                                    'Subgroup_25',
                                                                                    'Dyella',
                                                                                    'Rhodobium',
                                                                                    'Sandaracinus'))

#qPCR_T1_D60_processed$Compartment <- factor(qPCR_T1_D60_processed$Compartment, levels = c("BS", "RH", "RP"))

#Convert variables on the x-axis to factors with desired order
qPCR_T1_D60_processed$Genus_id <- factor(qPCR_T1_D60_processed$Genus_id, levels = rev(sort(unique(qPCR_T1_D60_processed$Genus_id))))

heatmap_qPCR_T1_D60<-ggplot(data=qPCR_T1_D60_processed,mapping=aes(x=Tests, y=Genus_id, fill=Result))
heatmap_qPCR_T1_D60 + 
  geom_tile(colour="darkgrey",size=0.25) +
  facet_grid(~Tests, scales="free", space="free") +
  theme(legend.position = "none") +
  xlab(label="") +
  scale_fill_gradient(name="Index",low="white", high="black", na.value="lightgrey") +
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  theme(axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size=16, colour="black", face = "italic")) +
  #theme(legend.title = element_text(size=16), legend.text = element_text(size = 16, colour = "black")) +
  theme(legend.title = element_blank()) +
  theme(strip.text.x = element_text(size = 16, colour = "black", angle = 90, hjust=0)) +
  theme(strip.text.y = element_text(size = 16, colour = "black", angle = 0)) +
  theme(strip.placement = "top", plot.title = element_text(hjust=0.5), axis.title.y = element_blank(), strip.background=element_rect(fill="#EEEEEE", color="#424242")) +
  ggtitle(label="")

ggsave("heatmap_qPCR_T1_D60.png", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Heatmap_jki_seq5/", width = 20, height = 14, units = "cm",dpi = 300)



### qPCR_T2_D60
#load data 
qPCR_T2_D60<-read_excel("~/Documents/R_analysis/jki_seq5/data_jki_seq5/qPCR_bacteria_T2_D60.xlsx")
qPCR_T2_D60

###Preparing tidyr format table (variables in columns, observation in rows, values in cells), collecting a set of column names and place them into a single value column
qPCR_T2_D60_processed<-gather(data=qPCR_T2_D60,key = Tests, value = Result,amoA:nxrAB)
head(qPCR_T2_D60_processed)

###Convert data to square root 
#vt_long$Sqrt.abundance<-sqrt(vt.long$Absolute_Abundance)

#####Constructing the ggplot heatmap

### Heatmap GO_96  
#Convert variables on the x-axis to factors with desired order
qPCR_T2_D60_processed$Tests <- factor(qPCR_T2_D60_processed$Tests, levels = c("amoA", "nifH", "nirS", "nosZ", "napA", "narG", "nasA", "nirB", "nirK", "norB", "nrfA", "nxrAB"))
qPCR_T2_D60_processed$Genus_id <- factor(qPCR_T2_D60_processed$Genus_id, levels = c('Verrucosispora',
                                                                                    'Caulobacter',
                                                                                    'Gemmatimonas',
                                                                                    'Streptomyces',
                                                                                    'Subgroup_25',
                                                                                    'Devosia',
                                                                                    'Fimbriimonadaceae',
                                                                                    "Devosiaceae"))

#qPCR_T2_D60_processed$Compartment <- factor(qPCR_T2_D60_processed$Compartment, levels = c("BS", "RH", "RP"))

#Convert variables on the x-axis to factors with desired order
qPCR_T2_D60_processed$Genus_id <- factor(qPCR_T2_D60_processed$Genus_id, levels = rev(sort(unique(qPCR_T2_D60_processed$Genus_id))))

heatmap_qPCR_T2_D60<-ggplot(data=qPCR_T2_D60_processed,mapping=aes(x=Tests, y=Genus_id, fill=Result))
heatmap_qPCR_T2_D60 + 
  geom_tile(colour="darkgrey",size=0.25) +
  facet_grid(~Tests, scales="free", space="free") +
  theme(legend.position = "none") +
  xlab(label="") +
  scale_fill_gradient(name="Index",low="white", high="black", na.value="lightgrey") +
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  theme(axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size=16, colour="black", face = "italic")) +
  #theme(legend.title = element_text(size=16), legend.text = element_text(size = 16, colour = "black")) +
  theme(legend.title = element_blank()) +
  theme(strip.text.x = element_text(size = 16, colour = "black", angle = 90, hjust=0)) +
  theme(strip.text.y = element_text(size = 16, colour = "black", angle = 0)) +
  theme(strip.placement = "top", plot.title = element_text(hjust=0.5), axis.title.y = element_blank(), strip.background=element_rect(fill="#EEEEEE", color="#424242")) +
  ggtitle(label="")

ggsave("heatmap_qPCR_T2_D60.png", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Heatmap_jki_seq5/", width = 18, height = 9, units = "cm",dpi = 300)

# The end! : )
