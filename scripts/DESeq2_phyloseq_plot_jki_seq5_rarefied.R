##Define taxa deferentially abundant using DESeq2 on phyloseq pipeline

#Loading package
library("phyloseq")
library("ggplot2")
library("DESeq2")
library("RColorBrewer")


# colors <- list(
    Phylum =  c("Acidobacteriota" = "#ffc331",
                "Actinobacteriota" = "#ff8b60",
                "Bacteroidota"  = "#Bff4be",
                "Campilobacterota" = "#800080",
                "Chloroflexi" = "#2b8f22",
                "Cyanobacteria" = "#185113",
                "Deferribacterota" = "#d4ff31",
                "Desulfobacterota" = "#964B00",
                "Firmicutes" = "#4B69F5",
                "Gemmatimonadota" = "#363781",
                "Latescibacterota" = "#4a2500",
                "Myxococcota" = "#FFC0CB",
                "NB1-j" = "#b5814c",
                "Nitrospirota" = "#1361cf",
                "Patescibacteria" = "#4b0096",
                "Planctomycetota" = "#9fc5e8",
                "Proteobacteria" = "#ca250f",
                "Verrucomicrobiota" = "#ced2ce",
                "Zixibacteria" = "#e89fc5")

    
### T1
    ## Control vs. Pex = only 1
    ## Control vs. SynCom = only 3
    ## Control vs. SynCom_Pex = only 1
    
### T2  
    ## Control vs. Pex = only 2
    ## Control vs. SynCom =  23
    ## Control vs. SynCom_Pex = only 15
    
### T2 
    ## Control vs. SynCom =  23
    
##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_T2_Control_SynCom_treat$log2FoldChange, sigtab_T2_Control_SynCom_treat$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_T2_Control_SynCom_treat$Phylum = factor(as.character(sigtab_T2_Control_SynCom_treat$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_T2_Control_SynCom_treat$log2FoldChange, sigtab_T2_Control_SynCom_treat$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_T2_Control_SynCom_treat$Annotation = factor(as.character(sigtab_T2_Control_SynCom_treat$Annotation), levels=names(x))

# Define color palette
color_phylum_T2 <- c("Actinobacteriota" = "#ff8b60",
                     "Firmicutes" = "#4B69F5",
                     "Proteobacteria" = "#ca250f",
                     "Verrucomicrobiota" = "#ced2ce")

plot_sigtab_T2_Control_SynCom_treat <- ggplot(sigtab_T2_Control_SynCom_treat, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=4) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5, size = 18, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 20, color = "black")) +
  theme(axis.text.y = element_text(size = 24, color = "black")) +
  theme(legend.title = element_text(size = 24, color = "black")) +
  theme(legend.text = element_text(size = 24, color = "black")) +
  #scale_color_manual(values = color_phylum_T2) +
  scale_color_manual(values = color_phylum_T2, breaks = c('Actinobacteriota',
                                                          'Firmicutes',
                                                          'Proteobacteria',
                                                          'Verrucomicrobiota')) +
  scale_x_continuous(limits = c(-8, 8), breaks = seq(-8, 8, by = 4), label = c("-8", "-4", "0", "4", "8")) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_T2_Control_SynCom_treat

ggsave("plot_sigtab_T2_Control_SynCom_treat.png", path = "~/Documents/R_analysis/igz_seq2/output_igz_seq2/DESeq2_plot_igz_seq2/", width = 25, height = 15, units = "cm", dpi = 300, device = "png")

## The end! : )

