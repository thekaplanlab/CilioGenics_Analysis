

source("functions.R")
#source("generate_objects.R")
library(stringr)
library(cluster)
library(readxl)
library(data.table)
library(dplyr)
library(geneName)

pwd <- "/home/kaplanlab/Desktop/Desktop_Items/MUSTAFA_PIR/BLAST_Cluster/Homo_sapiens/blast_folder2/results/"

newscores_50 <- filtering(pwd, 0.001, 50)

gower_dist_50 <- daisy(newscores_50[, 2:73], metric = "gower")
cluster_50 <- hclust(gower_dist_50, method = "ward.D2")

tree_50 <- cutree(cluster_50, k = 40)
aa_50 <- cbind(newscores_50[, c(1, 2)], tree_50) # Genes with clusters

aa_50_temp <- aa_50
colnames(aa_50_temp)[3] <- "cluster_num"

a <- as.data.frame(cbind(newscores_50, tree_50))
a[, 2:73] <- lapply(a[, 2:73], function(x) as.numeric(as.character(x)))

# dfx <- data.frame(cluster_num = 0, score = 0, f1 = 0)
# 
# for (i in 1:40) {
#   allg <- rbind(ciliaryGenes, negativeTrain)
#   as <- aa_50[aa_50$tree_50 == i, ]
#   allg$pred <- ifelse(allg$Gene_name %in% as$Gene_name, 1, 0)
#   allg$pred <- as.factor(allg$pred)
#   allg$is.ciliary <- as.factor(allg$is.ciliary)
# 
#   cnfs <- confusionMatrix(allg$pred, allg$is.ciliary, positive = "1", mode = "prec_recall")
#   f1score <- cnfs$byClass[["Precision"]]
# 
#   b <- a[a$tree_50 == i, 2:45]
#   c <- a[a$tree_50 == i, 46:73]
#   score1 <- (0.28 * sum(b) - 0.88 * sum(c)) / length(b[[1]])
#   dfx1 <- data.frame(cluster_num = i, score = score1, f1 = f1score)
#   dfx <- rbind(dfx, dfx1)
# }

## *******************

dfx<-aa_50
dfx$score<-ifelse(dfx$tree_50 %in% c(31,37), 1, 0)
cluster_scores<-left_join(hgnc_names[,2], dfx[,c(1,4)], by = c("Approved symbol" = "Gene_name"))
cluster_scores$score[is.na(cluster_scores$score)]<-0

## *******************

dfx$total <- dfx$score * dfx$f1

## ************************************************
dfx<-aa_50
colnames(dfx)[3]<-"cluster_num"
scores<-data.frame(cluster_num = c(31,37,8,2,17,39,16,3,30,10,23,24,1,4:7,9,11:15,18:22,25:29,32:36,38,40),
                   score2 = c(1,0.9,0.6,0.6,0.5,0.5,0.4,0.3,0.2,0.1,0.1,0.1,rep(0,28)))

scores<-data.frame(cluster_num = c(31,37,1:30,32:36,38:40),
                   score2 = c(1,1,rep(0,38)))

dfx<-left_join(dfx, scores, by = "cluster_num")
dfx$phylogenetic_scores<-normalization(dfx$score2)
cluster_scores<-left_join(hgnc_names[,2], dfx[,c(1,5)], by = c("Approved symbol" = "Gene_name"))
cluster_scores$phylogenetic_scores[is.na(cluster_scores$phylogenetic_scores)]<-0
## ************************************************

cluster_motif<-cluster_scores[cluster_scores$`Approved symbol` %in% pan_motifx$Gene_name,]
cluster_motif_top<-cluster_motif[cluster_motif$phylogenetic_scores > 0.7,]
cluster_top<-cluster_scores[cluster_scores$phylogenetic_scores > 0.7,]

motif_half<-motif_scores2[motif_scores2$motif_score > 0,]
cluster_motif_half<-cluster_scores[cluster_scores$`Approved symbol` %in% motif_half$Gene_name,]
cluster_motif_halftop<-cluster_motif_half[cluster_motif_half$phylogenetic_scores > 0.7,]

sc_half<-single_cell_scores[single_cell_scores$single_cell_score >= 0.4,]
sc_motif_half<-sc_half[sc_half$Gene_name %in% motif_half$Gene_name,]





#cluster_reg_scores <- fread("files/lr_scores_clusters.txt")

cluster_scores <- left_join(aa_50_temp, dfx[, c(1, 5)], by = "cluster_num")
cluster_scores$norm_score <- ifelse(cluster_scores$total >= 0, cluster_scores$total, 0)
cluster_scores$norm_score[is.na(cluster_scores$norm_score)] <- 0
cluster_scores$norm_score <- normalization0(cluster_scores$norm_score)
cluster_scores <- unique(cluster_scores[, c(1, 4)])
colnames(cluster_scores)[2] <- "phylogenetic_scores"

write.table(cluster_scores, "cluster_scores.txt", sep = "\t", row.names = FALSE, quote = FALSE)

write.table(cluster_scores, "cluster_scores_only31and37.txt", sep = "\t", row.names = FALSE, quote = FALSE)

