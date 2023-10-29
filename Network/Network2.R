### Interaction network ----

library(GGally)
library(network)
library(readxl)
library(sna)
setwd("/Users/oktayismailkaplan/Desktop/Manuscripts/Ciliogenics/Figures/Network")

all_scores2<-read.table("all_scores2.txt", sep = "\t", header = TRUE,  fill = TRUE)
ciliaryGenes<-read.table("ciliaryGenes.txt", sep = "\t", header = TRUE)
negativeCiliary<-read.table("negativeCiliary.txt", sep = "\t", header = TRUE)
ciliarynet<-read.table("protein_interaction_file.txt", sep = "\t", header = TRUE)

################### First Genes ######################
genes <- read_xlsx("Table_S8.xlsx")
genes <- genes[-1,]
glist<- genes$CilioGenics_top_500

#################################################

genes <- read_xlsx("Table_S9.xlsx")
genes <- genes[-1,]
glist<- genes$`Gold Standard Gene List`

#################################################

genes <- read_xlsx("Table_S10.xlsx")
genes <- genes[-1,-2]
glist<- genes$`Negative Gene List`



# **************************************************************
###### 1 # this will draw all edges between selected nodes
ciliarynet2<-ciliarynet[ciliarynet$Interactor_A %in% glist | ciliarynet$Interactor_B %in% glist,]
net<-unique(c(ciliarynet2$Interactor_A, ciliarynet2$Interactor_B))
netnet<-ciliarynet[ciliarynet$Interactor_A %in% net & ciliarynet$Interactor_B %in% net,]
ciliarynet2<-unique(rbind(ciliarynet2, netnet))

###### 2 # this will only draw edges between selected genes and their interactome
ciliarynet2<-ciliarynet[ciliarynet$Interactor_A %in% glist | ciliarynet$Interactor_B %in% glist,]

# **************************************************************

ciliarynet2$is.ciliary<-ifelse(ciliarynet2$Interactor_B %in% ciliaryGenes$Gene_name, "Ciliary", 
                               ifelse(ciliarynet2$Interactor_B %in% negativeCiliary$Gene_name, "Non-ciliary", "Unknown"))

cilsimpnet<-unique(ciliarynet2[,1:2])
cilsimpnet<-cilsimpnet[cilsimpnet$Interactor_A != cilsimpnet$Interactor_B,]

nw <- network(cilsimpnet, directed = TRUE, matrix.type = "edgelist", loops = TRUE)

candlist<-all_scores2$Gene_name[all_scores2$total_score_sums >= 1.2]

vertexnames<-network.vertex.names(nw)
nw %v% "is.ciliary" = ifelse(vertexnames %in% ciliaryGenes$Gene_name, "Ciliary", 
                             ifelse(vertexnames %in% negativeCiliary$Gene_name, "Non-ciliary",
                                    ifelse(vertexnames %in% candlist, "Candidate", "Unknown")))

set.seed(2)
ggnet2(nw, color = "is.ciliary", node.size = 2, label = glist, label.size = 4,
       palette = c("Ciliary" = "palegreen", "Non-ciliary" = "tomato", "Unknown" = "steelblue", "Candidate" = "yellow2"),
       mode = "kamadakawai")
