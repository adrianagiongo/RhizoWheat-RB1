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

    

##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_T1$log2FoldChange, sigtab_T1$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_T1$Phylum = factor(as.character(sigtab_T1$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_T1$log2FoldChange, sigtab_T1$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_T1$Annotation = factor(as.character(sigtab_T1$Annotation), levels=names(x))

# Define color palette
color_phylum_T1 <- c("Acidobacteriota" = "#ffc331",
                     "Actinobacteriota" = "#ff8b60",
                     "Bacteroidota"  = "#Bff4be",
                     "Firmicutes" = "#4B69F5",
                     "Gemmatimonadota" = "#363781",
                     "Nitrospirota" = "#1361cf",
                     "Proteobacteria" = "#ca250f",
                     "Verrucomicrobiota" = "#ced2ce")

plot_sigtab_T1 <- ggplot(sigtab_T1, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=4) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5, size = 18, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 20, color = "black")) +
  theme(axis.text.y = element_text(size = 24, color = "black")) +
  theme(legend.title = element_text(size = 24, color = "black")) +
  theme(legend.text = element_text(size = 24, color = "black")) +
  #scale_color_manual(values = color_phylum_T2) +
  scale_color_manual(values = color_phylum_T2, breaks = c('Bacteroidota',
                                                          'Firmicutes',
                                                          'Proteobacteria',
                                                          'Verrucomicrobiota')) +
  scale_x_continuous(limits = c(-30, 10), breaks = seq(-30, 10, by = 10), label = c("-30", "-20", "-10", "0", "10")) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_T1

ggsave("plot_sigtab_T1.png", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/DESeq2_plot_jki_seq5/", width = 25, height = 15, units = "cm", dpi = 300, device = "png")


##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_T2$log2FoldChange, sigtab_T2$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_T2$Phylum = factor(as.character(sigtab_T2$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_T2$log2FoldChange, sigtab_T2$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_T2$Annotation = factor(as.character(sigtab_T2$Annotation), levels=names(x))

# Define color palette
color_phylum_T2 <- c("Acidobacteriota" = "#ffc331",
                     "Actinobacteriota" = "#ff8b60",
                     "Bacteroidota"  = "#Bff4be",
                     "Firmicutes" = "#4B69F5",
                     "Gemmatimonadota" = "#363781",
                     "Nitrospirota" = "#1361cf",
                     "Proteobacteria" = "#ca250f",
                     "Verrucomicrobiota" = "#ced2ce")

plot_sigtab_T2 <- ggplot(sigtab_T2, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=4) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5, size = 18, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 20, color = "black")) +
  theme(axis.text.y = element_text(size = 24, color = "black")) +
  theme(legend.title = element_text(size = 24, color = "black")) +
  theme(legend.text = element_text(size = 24, color = "black")) +
  #scale_color_manual(values = color_phylum_T2) +
  scale_color_manual(values = color_phylum_T2, breaks = c('Acidobacteriota',
                                                          'Actinobacteriota',
                                                          'Bacteroidota',
                                                          'Firmicutes',
                                                          'Gemmatimonadota',
                                                          'Nitrospira',
                                                          'Proteobacteria',
                                                          'Verrucomicrobiota')) +
  scale_x_continuous(limits = c(-10, 30), breaks = seq(-10, 30, by = 10), label = c("-10", "0", "10", "20", "30")) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_T2

ggsave("plot_sigtab_T2.png", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/DESeq2_plot_jki_seq5/", width = 25, height = 16, units = "cm", dpi = 300, device = "png")


## The end! : )

