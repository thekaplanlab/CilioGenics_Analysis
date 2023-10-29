
## required libraries ##

library(readxl)
library(dplyr)
library(stringr)

## set working directory ##

setwd("C:\\Users\\Hasan Can Demirci\\Desktop\\CilioGenics Revision\\count_data_motile")

## uploud file ##

all <- read.csv2("all_genes_scores.csv")

## select and manipulate data ##

selected_columns <- all %>%
  select(Reyfman, Carraro, Habermann, Murthy, Trachea) %>%
  mutate(Counts = rowSums(!is.na(.)))

selected_columns <- data.frame(selected_columns, all$Gene_name)

selected_columns <- selected_columns[c(7,1,2,3,4,5,6)]

names(selected_columns)[1] <- "Gene_name"

## write column for having gene ##

selected_columns$Studies <- apply(selected_columns[, 2:6], 1, function(row) {
  non_na_columns <- names(row[row > 0])
  if (length(non_na_columns) > 0) {
    paste(non_na_columns, collapse = "   ")
  } else {
    ""
  }
})

## NA to empty ##

selected_columns$Studies <- str_replace_all(selected_columns$Studies, "NA", "    ")

## select and filter ##

selected_columns <- selected_columns %>%
  select(1,7,8) %>%
  filter(Counts >0)

## export data to working directory ##

write.csv(selected_columns,"count_data_motile.csv")






