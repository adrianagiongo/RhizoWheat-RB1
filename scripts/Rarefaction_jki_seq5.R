#Create rarefied data from phyloseq object using vegan package or microbiome package
# https://rdrr.io/cran/vegan/man/rarefy.html#google_vignette

#Loading package
library("phyloseq")
library("vegan")
library("microbiome")
library("ggpubr")
library("ggplot2")
library("tidyr")

##Creating rarefaction curve of non-rarefied samples
#transpose S4 otu-table
mat_psO_jki_seq5 <- t(otu_table(psO_jki_seq5))

#Transform S4 objet to matrix (a warning message is normal)
class(mat_psO_jki_seq5) <- "matrix"

#Create curve using rarecurve()
system.time(rarecurve(mat_psO_jki_seq5, step = 1000, lwd=2, ylab="OTU", col = "black", label = FALSE))

##Rarefying samples to the minimum number of reads among samples
#Using rarefy_even_depth() rarefy to the lower number of total sequences in a sample
psO_jki_seq5_rarefied<-rarefy_even_depth(psO_jki_seq5, rngseed=2022, sample.size = min(sample_sums(psO_jki_seq5)),trimOTUs=TRUE)
psO_jki_seq5_rarefied
sample_sums(psO_jki_seq5_rarefied)

##Creating rarefaction curve of non-rarefied samples
#transpose S4 otu-table
mat_psO_jki_seq5_rarefied <- t(otu_table(psO_jki_seq5_rarefied))

#Transform S4 objet to matrix
class(mat_psO_jki_seq5_rarefied) <- "matrix"

#Create curve using rarecurve()
system.time(rarecurve(mat_psO_jki_seq5_rarefied, step = 1000, lwd=2, ylab="OTU", col = "blue", label = FALSE))

#Prepare file with meta data using meta function
psO_jki_seq5_rarefied_rare.meta <- meta(psO_jki_seq5_rarefied)
head(psO_jki_seq5_rarefied_rare.meta)
psO_jki_seq5_rarefied_rare.meta.rot_mic <- subset(psO_jki_seq5_rarefied_rare.meta, select = -c(Sample_name, Location, Rotation, Stage, Microhabitat, Depth, Replicate, Rot_Stage, Rot_Stage_Mic))
head(psO_jki_seq5_rarefied_rare.meta.rot_mic)

#Create and transpose matrix to make x axis for samples
df_psO_jki_seq5_rarefied <- data.frame(otu_table(psO_jki_seq5_rarefied))
#write.csv(df_psO_jki_seq5_rarefied, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/df_psO_jki_seq5_rarefied.csv")

df_psO_jki_seq5_rarefied_t <- t(df_psO_jki_seq5_rarefied)

#Combine matrix in one
ASVs_rarefied_ed<-cbind(psO_jki_seq5_rarefied_rare.meta.rot_mic, df_psO_jki_seq5_rarefied_t)
ASVs_rarefied_ed[1:36,1]

##check for the distribution of any ASV counts in the sample group using "hist" function to plot, 
#"shapiro.test" function to test the Null hypothesis, and "qqnorm" function to qq plot
# If alpha result is lower than the alpha value chosen (p<0.05) 
#then the null hypotesis (the population is normally distributed) is rejected
hist(ASVs_rarefied_ed$sp1)
shapiro.test(ASVs_rarefied_ed$sp1)
qqnorm(ASVs_rarefied_ed$sp1)

#Define colours
colour = rep(NA, length=length(ASVs_rarefied_ed[,1]))
colour[which(ASVs_rarefied_ed$Rot_Mic=="W1_RA")] = "#004d26"
colour[which(ASVs_rarefied_ed$Rot_Mic=="W1_RH")] = "#00b33c"
colour[which(ASVs_rarefied_ed$Rot_Mic=="W2_RA")] = "#b22222"
colour[which(ASVs_rarefied_ed$Rot_Mic=="W2_RH")] = "#f4837c"
dim(ASVs_rarefied_ed)

#plot
tiff("~/Documents/R_analysis/jki_seq5/output_jki_seq5/Alpha_div_jki_seq5/rarecurve_jki_seq5_rarefied.tiff", units="cm", width=12, height=12, res=300)
rarecurve(ASVs_rarefied_ed [,2:15862], step=1000, label=FALSE, col=colour, main="", xlab = "Number of sequences", ylab = "ASVs")
legend(legend=c("W1_RA", "W1_RH", "W2_RA", "W2_RH"),
       "bottomright", 
       bty="n",
       col=c("#004d26", "#00b33c", "#b22222", "#f4837c"),
       pch=15,
       pt.cex=1.5,
       cex=0.75,
       ncol=3)
#arrows(x0=30464, y0=300, y1=1800, angle=90, length=0, lty = 2)
dev.off()


## The end! Have fun! Enjoy the day!  : )
