library(data.table)
library(readxl)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(colorspace)

#cilia_TFs <- as.table(rbind(c(562, 270, 246, 878, 126, 433), c(127, 89, 99, 178, 50, 107)))
#dimnames(cilia_TFs) <- list(Genes = c("Candidate Genes", "Ciliary Genes"),
                                    #TZ = c("GLIS3","JAZF1", "RFX2", "RFX3", "MYB", "FOXJ1"))

#library
library(ggplot2)
library(viridis)
library(readxl)
setwd("/Users/oktay/Desktop/Manuscripts/Ciliogenics/Figures")
cilia_TFs <- read_excel("TZ.xlsx")

View(cilia_TFs)
cilia_TFs

x <- ggplot(cilia_TFs, aes(fill=Genes, y=Number, x=Gene_Name))+
  geom_bar(position="fill", stat="identity")+
  scale_fill_viridis(discrete = T) +
  ggtitle("Transcription Factors") +
  theme_classic()+
  scale_fill_brewer(palette = 1)

x + scale_x_discrete(limits = c("GLIS3", "JAZF1", "RFX2", "RFX3", "MYB","FOXJ1"))