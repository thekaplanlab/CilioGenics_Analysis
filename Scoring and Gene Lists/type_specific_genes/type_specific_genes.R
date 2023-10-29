
## required library ##

library(readxl)
library(dplyr)
library(stringr)

## set working directory ##

setwd("C:\\Users\\Hasan Can Demirci\\Desktop\\CilioGenics Revision\\type_specific_genes")

## upload the files ##

all <- read.csv2("all_genes_scores.csv")

## select motile columns ##

motilGenes <- all[c(1,3,4,5,6,8)]

## filter depends on scores ##

m_filtered_data <- motilGenes %>%
  filter(
    Reyfman >=0,
    Carraro >=0,
    Habermann >=0,
    Murthy >=0,
    Trachea >=0)

## export the last file to working directory ##

write.csv2(m_filtered_data, "motile_genes.csv")


## upload the files for primary cilia ##

Cengen <- read.csv2("cengen_all_scores_human.csv")
Cao <- read_xlsx("Cao(2017).xlsx")
combined_primary <- read_xlsx("combined_cao_cengen.xlsx")

## select the target column ##

Cengen <- Cengen %>%
  select(13) 

## make unique data values ##

Cao <- unique(Cao)

## Score the genes depends on existence ##

combined_primary$Cao <- ifelse(combined_primary$genes %in% Cao$`Cao(2017)`,1,0)

combined_primary$CenGen <- ifelse(combined_primary$genes %in% Cengen$human_symbol,1,0)

## make unique the data values ##

combined_primary <- unique(combined_primary)

## filter data based on existence score ##

combined_primary <- combined_primary%>%
  filter(Cao == 1 & CenGen == 1)

## export the last file to working directory ##

write.csv(combined_primary,"primary_genes.csv")


                   