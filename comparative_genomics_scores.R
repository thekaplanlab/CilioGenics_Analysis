
library(dplyr)
library(cluster)
library(geneName)

source("functions.R")

## Comparative genomics ----

pwd <- "././files/comparative_genomics"

newscores_50 <- filtering(pwd, 0.001, 50)

gower_dist_50 <- daisy(newscores_50[, 2:73], metric = "gower")
cluster_50 <- hclust(gower_dist_50, method = "ward.D2")

tree_50 <- cutree(cluster_50, k = 40)
aa_50 <- cbind(newscores_50[, c(1, 2)], tree_50) # Genes with clusters

aa_50_temp <- aa_50
colnames(aa_50_temp)[3] <- "cluster_num"

a <- as.data.frame(cbind(newscores_50, tree_50))
a[, 2:73] <- lapply(a[, 2:73], function(x) as.numeric(as.character(x)))

dfx<-aa_50
colnames(dfx)[3]<-"cluster_num"
scores<-data.frame(cluster_num = c(31,37,1:30,32:36,38:40),
                   score2 = c(1,1,rep(0,38)))


dfx<-hgncConverter(dfx, "Gene_name")
dfx<-left_join(dfx, scores, by = "cluster_num")
dfx$phylogenetic_scores<-normalization(dfx$score2)
cluster_scores<-left_join(hgnc_names[,2], dfx[,c(1,5)], by = c("Approved symbol" = "Gene_name"))
cluster_scores$phylogenetic_scores[is.na(cluster_scores$phylogenetic_scores)]<-0
cluster_scores<-cluster_scores %>%
  group_by(`Approved symbol`) %>%
  summarise_all(max)
colnames(cluster_scores)[1]<-"Gene_name"
write.table(cluster_scores, "cluster_scores.txt", sep = "\t", row.names = FALSE, quote = FALSE)


