#ANCOM-BC
# https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html
#ANCOM-BC2
# https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html

library(ANCOMBC)
library(microbiome)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(DT)
options(DT.options = list(
  initComplete = JS("function(settings, json) {",
                    "$(this.api().table().header()).css({'background-color': 
  '#000', 'color': '#fff'});","}")))
library("phyloseq")
library("readxl")       
library("tibble")       
library("tidyverse")   
library("vegan")
library("ggpubr")
library("RColorBrewer")

####ANCOM-BC implementation
#we consider the following covariates:
#Continuous covariates: none, could be for example: “age”
#Categorical covariates: “site”, “rootstock”
#The group variable of interest: “site”
#Three groups: “EH”, “HG”, “RU”
#The reference group: “EH”

#by default, levels of a categorical variable in R are sorted alphabetically. In this case, the reference level for `site` will be 
# `EH`. To manually change the reference level, for instance, setting `HG` as the reference level, use:
#df$site = factor(df$site, levels = c("HG", "EH", "RU"))

#run anbombc2 (importing data in phyloseq format)
set.seed(2608)
output = ancombc2(data = psO_jki_seq5_rarefied_T1_filt_annotation, assay_name = "counts", tax_level = "Family",
                  fix_formula = "Rotation + Microhabitat", rand_formula = NULL,
                  p_adj_method = "hochberg", 
                  group = "Rotation", 
                  alpha = 0.05, n_cl = 2, 
                  struc_zero = TRUE,
                  global = TRUE, pairwise = TRUE, 
                  dunnet = TRUE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = NULL, 
                  mdfdr_control = list(fwer_ctrl_method = "hochberg", B = 100), 
                  trend_control = NULL)

#detection of structural zeros (taxon presence/absence)
tab_zero = output$zero_ind
tab_zero %>%
  datatable(caption = "The detection of structural zeros")
View(tab_zero)
write.table(tab_zero, file = "tab_zero.csv", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

res_prim = output$res
View(res_prim)
write.table(res_prim, file = "WP3_ITS_ANCOM-BC2_annot_raref_res_prim.csv", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)


#ANCOMBC primary result
#Result from the ANCOM-BC log-linear model to determine taxa that are differentially abundant according to the covariate of interest. 
#It contains: 
#1) log fold changes; 
#2) standard errors; 
#3) test statistics; 
#4) p-values; 
#5) adjusted p-values; 
#6) indicators whether the taxon is differentially abundant (TRUE) or not (FALSE).



##ANCOM-BC2 global test and heat map of taxa with differential abundance across the “RU”, “HG”, and “EH” categories.
#In the subsequent heatmap, each cell represents a log fold-change (in natural log) value. Taxa marked in green have successfully passed the sensitivity analysis for pseudo-count addition.
res_global = output$res_global
View(res_global)
write.table(res_global, file = "WP3_ITS_ANCOM-BC2_annot_raref_res_global.csv", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
df_site = res_prim %>%
  dplyr::select(taxon, contains("site")) 
df_fig_global = df_site %>%
  dplyr::left_join(res_global %>%
                     dplyr::transmute(taxon, 
                                      diff_site = diff_abn, 
                                      passed_ss = passed_ss)) %>%
  dplyr::filter(diff_site == 1) %>%
  dplyr::mutate(lfc_HG = lfc_siteHG,
                lfc_RU = lfc_siteRU,
                color = ifelse(passed_ss == 1, "aquamarine3", "black")) %>%
  dplyr::transmute(taxon,
                   `HG - EH` = round(lfc_HG, 2),
                   `RU - EH` = round(lfc_RU, 2), 
                   color = color) %>%
  tidyr::pivot_longer(cols = `HG - EH`:`RU - EH`, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)

df_fig_global$group = factor(df_fig_global$group, 
                             levels = c("HG - EH",
                                        "RU - EH"))

lo = floor(min(df_fig_global$value))
up = ceiling(max(df_fig_global$value))
mid = (lo + up)/2
fig_global = df_fig_global %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Log fold changes for globally significant taxa") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(color = df_fig_global %>%
                                     dplyr::distinct(taxon, color) %>%
                                     .$color))
fig_global


