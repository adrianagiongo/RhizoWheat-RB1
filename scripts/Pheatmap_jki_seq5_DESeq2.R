#load packages
library("pheatmap")
library("dendextend")
library("vegan")
library("RColorBrewer")

## T1_D30
#load table
Pheatmap_jki_seq5_rarefied_T1_D30_data<-read.csv("~/Documents/R_analysis/jki_seq5/data_jki_seq5/Pheatmap/Pheatmap_jki_seq5_rarefied_T1_D30_filt_annotation_rel_DESeq2.csv", sep= ",", header = TRUE, row.names = 1)
dim(Pheatmap_jki_seq5_rarefied_T1_D30_data)
head(Pheatmap_jki_seq5_rarefied_T1_D30_data[,1:10])

#transpose matrix to make x axis for samples
Pheatmap_jki_seq5_rarefied_T1_D30_data_t <- t(Pheatmap_jki_seq5_rarefied_T1_D30_data)

#load metadata
#for columns
Pheatmap_jki_seq5_rarefied_col_meta<-read.csv("~/Documents/R_analysis/jki_seq5/data_jki_seq5/Pheatmap/Pheatmap_jki_seq5_rarefied_T1_D30_filt_annotation_rel_DESeq2_col.csv", sep= ",", header = TRUE, row.names = 1)
dim(Pheatmap_jki_seq5_rarefied_col_meta)
str(Pheatmap_jki_seq5_rarefied_col_meta)

#define colors
colors <- list(
  Rotation = c("W1" = "#d7d3a9", "W2" = "#74a553"),
  Microhabitat = c("RA" = "#b05644", "RH" = "#d9b967"))
 
#create heatmap using pheatmap and save
#Set colors
paletteLength <- 100
jki_seq5_colors <- colorRampPalette(c("white", "black"))(paletteLength)

#RA_differential_data
plot_Pheatmap_jki_seq5_rarefied_T1_D30_data_t<-pheatmap(Pheatmap_jki_seq5_rarefied_T1_D30_data_t, color = jki_seq5_colors, 
                                               border_color = "black", cluster_cols = TRUE, cluster_rows = FALSE,
                                               cellwidth = 8, cellheight = 10, annotation_col = Pheatmap_jki_seq5_rarefied_col_meta,
                                               annotation_colors = colors, angle_col = "90", cex=1.0, 
                                               fontsize_col = 8, fontsize_row = 10)

save_plot_Pheatmap_jki_seq5_rarefied_T1_D30_data_t <- function(x, filename, width=1900, height= 1400, res = 300) {
  tiff(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_plot_Pheatmap_jki_seq5_rarefied_T1_D30_data_t(plot_Pheatmap_jki_seq5_rarefied_T1_D30_data_t, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Heatmap_jki_seq5/save_plot_Pheatmap_jki_seq5_rarefied_T1_D30_data_t.png")

########################################
## T2_D30
#load table
Pheatmap_jki_seq5_rarefied_T2_D30_data<-read.csv("~/Documents/R_analysis/jki_seq5/data_jki_seq5/Pheatmap/Pheatmap_jki_seq5_rarefied_T2_D30_filt_annotation_rel_DESeq2.csv", sep= ",", header = TRUE, row.names = 1)
dim(Pheatmap_jki_seq5_rarefied_T2_D30_data)
head(Pheatmap_jki_seq5_rarefied_T2_D30_data[,1:5])

#transpose matrix to make x axis for samples
Pheatmap_jki_seq5_rarefied_T2_D30_data_t <- t(Pheatmap_jki_seq5_rarefied_T2_D30_data)

#load metadata
#for columns
Pheatmap_jki_seq5_rarefied_col_meta<-read.csv("~/Documents/R_analysis/jki_seq5/data_jki_seq5/Pheatmap/Pheatmap_jki_seq5_rarefied_T2_D30_filt_annotation_rel_DESeq2_col.csv", sep= ",", header = TRUE, row.names = 1)
dim(Pheatmap_jki_seq5_rarefied_col_meta)
str(Pheatmap_jki_seq5_rarefied_col_meta)

#define colors
colors <- list(
  Rotation = c("W1" = "#d7d3a9", "W2" = "#74a553"),
  Microhabitat = c("RH" = "#d9b967"))

#create heatmap using pheatmap and save
#Set colors
paletteLength <- 100
jki_seq5_colors <- colorRampPalette(c("white", "#0CAFA9"))(paletteLength)

#RA_differential_data
plot_Pheatmap_jki_seq5_rarefied_T2_D30_data_t<-pheatmap(Pheatmap_jki_seq5_rarefied_T2_D30_data_t, color = jki_seq5_colors, 
                                                    border_color = "black", cluster_cols = TRUE, cluster_rows = FALSE,
                                                    cellwidth = 8, cellheight = 10, annotation_col = Pheatmap_jki_seq5_rarefied_col_meta,
                                                    annotation_colors = colors, angle_col = "90", cex=1.0, 
                                                    fontsize_col = 8, fontsize_row = 10)

save_plot_Pheatmap_jki_seq5_rarefied_T2_D30_data_t <- function(x, filename, width=1900, height= 1400, res = 300) {
  tiff(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_plot_Pheatmap_jki_seq5_rarefied_T2_D30_data_t(plot_Pheatmap_jki_seq5_rarefied_T2_D30_data_t, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Heatmap_jki_seq5/save_plot_Pheatmap_jki_seq5_rarefied_T2_D30_data_t.tiff")


## T1_D60
#load table
Pheatmap_jki_seq5_rarefied_T1_D60_data<-read.csv("~/Documents/R_analysis/jki_seq5/data_jki_seq5/Pheatmap/Pheatmap_jki_seq5_rarefied_T1_D60_filt_annotation_rel_DESeq2.csv", sep= ",", header = TRUE, row.names = 1)
dim(Pheatmap_jki_seq5_rarefied_T1_D60_data)
head(Pheatmap_jki_seq5_rarefied_T1_D60_data[,1:10])

#transpose matrix to make x axis for samples
Pheatmap_jki_seq5_rarefied_T1_D60_data_t <- t(Pheatmap_jki_seq5_rarefied_T1_D60_data)

#load metadata
#for columns
Pheatmap_jki_seq5_rarefied_col_meta<-read.csv("~/Documents/R_analysis/jki_seq5/data_jki_seq5/Pheatmap/Pheatmap_jki_seq5_rarefied_T1_D60_filt_annotation_rel_DESeq2_col.csv", sep= ",", header = TRUE, row.names = 1)
dim(Pheatmap_jki_seq5_rarefied_col_meta)
str(Pheatmap_jki_seq5_rarefied_col_meta)

#define colors
colors <- list(
  Rotation = c("W1" = "#d7d3a9", "W2" = "#74a553"),
  Microhabitat = c("RA" = "#b05644", "RH" = "#d9b967"))

#create heatmap using pheatmap and save
#Set colors
paletteLength <- 100
jki_seq5_colors <- colorRampPalette(c("white", "#0CAFA9"))(paletteLength)

#RA_differential_data
plot_Pheatmap_jki_seq5_rarefied_T1_D60_data_t<-pheatmap(Pheatmap_jki_seq5_rarefied_T1_D60_data_t, color = jki_seq5_colors, 
                                                        border_color = "black", cluster_cols = TRUE, cluster_rows = FALSE,
                                                        cellwidth = 8, cellheight = 10, annotation_col = Pheatmap_jki_seq5_rarefied_col_meta,
                                                        annotation_colors = colors, angle_col = "90", cex=1.0, 
                                                        fontsize_col = 8, fontsize_row = 10)

save_plot_Pheatmap_jki_seq5_rarefied_T1_D60_data_t <- function(x, filename, width=1900, height= 1400, res = 300) {
  tiff(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_plot_Pheatmap_jki_seq5_rarefied_T1_D60_data_t(plot_Pheatmap_jki_seq5_rarefied_T1_D60_data_t, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Heatmap_jki_seq5/save_plot_Pheatmap_jki_seq5_rarefied_T1_D60_data_t.tiff")

########################################
## T2_D60
#load table
Pheatmap_jki_seq5_rarefied_T2_D60_data<-read.csv("~/Documents/R_analysis/jki_seq5/data_jki_seq5/Pheatmap/Pheatmap_jki_seq5_rarefied_T2_D60_filt_annotation_rel_DESeq2.csv", sep= ",", header = TRUE, row.names = 1)
dim(Pheatmap_jki_seq5_rarefied_T2_D60_data)
head(Pheatmap_jki_seq5_rarefied_T2_D60_data[,1:5])

#transpose matrix to make x axis for samples
Pheatmap_jki_seq5_rarefied_T2_D60_data_t <- t(Pheatmap_jki_seq5_rarefied_T2_D60_data)

#load metadata
#for columns
Pheatmap_jki_seq5_rarefied_col_meta<-read.csv("~/Documents/R_analysis/jki_seq5/data_jki_seq5/Pheatmap/Pheatmap_jki_seq5_rarefied_T2_D60_filt_annotation_rel_DESeq2_col.csv", sep= ",", header = TRUE, row.names = 1)
dim(Pheatmap_jki_seq5_rarefied_col_meta)
str(Pheatmap_jki_seq5_rarefied_col_meta)

#define colors
colors <- list(
  Rotation = c("W1" = "#d7d3a9", "W2" = "#74a553"),
  Microhabitat = c("RH" = "#d9b967"))

#create heatmap using pheatmap and save
#Set colors
paletteLength <- 100
jki_seq5_colors <- colorRampPalette(c("white", "#0CAFA9"))(paletteLength)

#RA_differential_data
plot_Pheatmap_jki_seq5_rarefied_T2_D60_data_t<-pheatmap(Pheatmap_jki_seq5_rarefied_T2_D60_data_t, color = jki_seq5_colors, 
                                                        border_color = "black", cluster_cols = TRUE, cluster_rows = FALSE,
                                                        cellwidth = 8, cellheight = 10, annotation_col = Pheatmap_jki_seq5_rarefied_col_meta,
                                                        annotation_colors = colors, angle_col = "90", cex=1.0, 
                                                        fontsize_col = 8, fontsize_row = 10)

save_plot_Pheatmap_jki_seq5_rarefied_T2_D60_data_t <- function(x, filename, width=1900, height= 1400, res = 300) {
  tiff(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_plot_Pheatmap_jki_seq5_rarefied_T2_D60_data_t(plot_Pheatmap_jki_seq5_rarefied_T2_D60_data_t, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Heatmap_jki_seq5/save_plot_Pheatmap_jki_seq5_rarefied_T2_D60_data_t.tiff")



### Enjoy the sun, have fun! : )
