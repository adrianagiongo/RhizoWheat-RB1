##Supervised cleaning process to clean the raw dataset

#load packages
library("phyloseq")
library("ggplot2")
library("dplyr")

####################################### 
#Show available ranks in the dataset  --> original dataset = 153 Go + 136 Ki = 289 samples
rank_names(jki_seq5_data)
sample_names(jki_seq5_data)
jki_seq5_data

##Create table with number of features for each phylum
##Remove respective artefacts on a specific taxonomic level (Keeping Archaea)
#Kingdom
table(tax_table(jki_seq5_data)[,"Kingdom"], exclude = NULL)
psO_jki_seq5 <- subset_taxa(jki_seq5_data, !is.na(Kingdom) & !Kingdom %in% c("Unassigned", "Eukaryota"))
table(tax_table(psO_jki_seq5)[,"Kingdom"], exclude = NULL)

#Phylum
table(tax_table(psO_jki_seq5)[,"Phylum"], exclude = NULL)
psO_jki_seq5 <- subset_taxa(psO_jki_seq5, !Phylum %in% c("Archaea")) #No NAs anymore (Archaea instead)
table(tax_table(psO_jki_seq5)[,"Phylum"], exclude = NULL)

#Order
table(tax_table(psO_jki_seq5)[,"Order"], exclude = NULL)
psO_jki_seq5 <- subset_taxa(psO_jki_seq5, !is.na(Order) & !Order %in% c("Chloroplast"))
table(tax_table(psO_jki_seq5)[,"Order"], exclude = NULL)

#Family
table(tax_table(psO_jki_seq5)[,"Family"], exclude = NULL)
psO_jki_seq5 <- subset_taxa(psO_jki_seq5, !is.na(Family) & !Family %in% c("Mitochondria"))
table(tax_table(psO_jki_seq5)[,"Family"], exclude = NULL)

#Genus
#table(tax_table(psO_jki_seq5)[,"Genus"], exclude = NULL)

#Annotation
#table(tax_table(psO_jki_seq5)[,"Annotation"], exclude = NULL)

##Observe psO after clean expurious taxa
jki_seq5_data
psO_jki_seq5

#Compute prevalence of each feature and store as data.frame
prevdf = apply(X = otu_table(psO_jki_seq5),
               MARGIN = ifelse(taxa_are_rows(psO_jki_seq5), yes = 1, no = 2),
               FUN = function (x) {sum(x>0)})

#Add taxonomy and total reads counts to the data frame
prevdf_seq5 = data.frame(Prevalence=prevdf, TotalAbundance = taxa_sums(psO_jki_seq5), tax_table(psO_jki_seq5), otu_table(psO_jki_seq5))

#Write data frame in csv format
write.csv(prevdf_seq5, "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/psO_jki_seq5_prevdf.csv")


## The end  : )
