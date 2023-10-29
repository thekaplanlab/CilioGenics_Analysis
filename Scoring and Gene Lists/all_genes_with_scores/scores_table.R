
## required libraries ##

library(readxl)
library(dplyr)

## set working directory ##

setwd("C:\\Users\\Hasan Can Demirci\\Desktop\\CilioGenics Revision\\all_genes_with_scores")

## uploading files ##

trachea <- read.csv2("trachea_ciliated_markers.csv",sep = ",")

Cengen <- read.csv2("cengen_all_scores_human.csv", sep = ",")
all <- read_xlsx("o_all_scores_sidebyside.xlsx")

## Select only sensory ones ##

Cengen <- Cengen %>%
  filter(cluster == "Sensory")

## new Cengen data ##

write.csv2(Cengen,"cengen_all_scores_human.csv")

## creating new columns for new analysis and matching their scores ##

all$Trachea <- ifelse(all$Gene_name %in% trachea$gene_name,
                      trachea$score[match(all$Gene_name ,trachea$gene_name)],
                      NA)

all$CenGen <- ifelse(all$Gene_name %in% Cengen$human_symbol,
                     Cengen$score[match(all$Gene_name ,Cengen$human_symbol)],
                     NA)

## exporting the last scores table to working directory ##

write.csv2(all,"all_genes_scores2.csv")
