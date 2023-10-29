
## required libraries ##

library(VennDiagram)
library(gridExtra)
library(readxl)

## set working directory ##

setwd("C:\\Users\\Hasan Can Demirci\\Desktop\\CilioGenics Revision\\gene_extraction_Venn")

## upload files ##

motile <- read.csv2("motile_genes.csv")
primary <- read.csv("primary_genes.csv")
gold <- read_xlsx("gold_genes.xlsx")
negative <- read.table("negativeCiliary.txt", header = T)

## select genes column ##

motile <- motile[1]
primary <- primary[2]
negative <- negative[1]

## other ones already unique, make unique negative ##

negative <- unique(negative)

## convert character format for Venn ##

motile <- as.character(motile$Gene_name)
primary <- as.character(primary$genes)
gold <- as.character(gold$Gene_name)
negative <- as.character(negative$Gene_name)

## get venn plot ##

venn <- venn.diagram(
  x = list(motile, primary,gold,negative),
  category.names = c("Motile", "Primary","Gold Standard","Negative Ciliary"),
  filename = NULL,
  output = TRUE,
  fill = c("red", "blue","green","orange")
)

grid.arrange(venn)

## Getting all intersections ##

library(gplots) #required libray 

## create list format for all data ##

listtt <-list(motile, primary,gold,negative)

## get plot by "TRUE" ##
ItemsList <- venn(listtt, show.plot = FALSE)

## check the dimensions of intersections ##

lengths(attributes(ItemsList)$intersections)

## check all members depends on intersections ##

attributes(ItemsList)$intersections

## convert intersection members to list format,then, matrix, after data_frame ##

last_Data <- list(attributes(ItemsList))

matrixx <- matrix(last_Data[[1]][["intersections"]])

matrix_data <- data.frame(matrixx)

## extract all intersection members step by step ##

motile <- data.frame(matrix_data[1,])
names(motile)[1]<-"motile_Genes"

primary <- data.frame(matrix_data[2,])
names(primary)[1]<-"primary_Genes"

gold <- data.frame(matrix_data[3,])
names(gold)[1]<-"gold_Genes"

negative <- data.frame(matrix_data[4,])
names(negative)[1]<-"negative_Genes"

motile_primary <- data.frame(matrix_data[5,])
names(motile_primary)[1]<-"motile_primary_Genes"

motile_gold <- data.frame(matrix_data[6,])
names(motile_gold)[1]<-"motile_gold_Genes"

motile_negative <- data.frame(matrix_data[7,])
names(motile_negative)[1]<-"motile_negative_Genes"

primary_gold <- data.frame(matrix_data[8,])
names(primary_gold)[1]<-"primary_gold_Genes"

primary_negative <- data.frame(matrix_data[9,])
names(primary_negative)[1]<-"primary_negative_Genes"

motile_gold_primary <- data.frame(matrix_data[10,])
names(motile_gold_primary)[1]<-"motile_gold_primary_Genes"

## combined all results ##

all_cilia_geness <-  bind_rows(motile, primary, gold, negative,
                               motile_primary, motile_gold, motile_negative,primary_gold,primary_negative,motile_gold_primary)

## extract Na and non-NA to order the values ##

all_cilia_geness <- lapply(all_cilia_geness, function(x) c(x[!is.na(x)], x[is.na(x)]))
all_cilia_geness <- as.data.frame(all_cilia_geness)


## export data to working directory ##

write.csv(all_cilia_geness,"primary_motile_gold.csv")




















