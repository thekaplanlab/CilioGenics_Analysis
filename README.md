**CilioGenics is an integrated method for predicting the ciliary genes**

To more accurately predict ciliary genes, CilioGenics combines several approaches, including as single-cell RNA sequencing, protein-protein interactions (PPIs), comparative genomics, transcription factor (TF)-network analysis, and text mining. 

Due to space limitations, the current repository contains programs for each strategy.


**Single cell RNA-seq analysis code**

**Load the following package**


``` Python
library(Seurat)
```

Dowload scRNA-seq file and read it 

``` Python
reyfman_dir<-dir("./sc_files/Reyfman")
reyfman_dir<-reyfman_dir[grepl("filtered", reyfman_dir)]
reyfman_dir_don<-reyfman_dir[grepl("Donor", reyfman_dir)]
reyfman_dir_dis<-reyfman_dir[!grepl("Donor", reyfman_dir)]
reyfman_dir_dis<-reyfman_dir_dis[!grepl("Cryobiopsy", reyfman_dir_dis)]
reyfman_all<-c(reyfman_dir_don, reyfman_dir_dis)
reyfmans<-list()
```

**Data visulaziation**
``` Python
ElbowPlot(reyfmans.all.combined, ndims = 30)
reyfmans.all.combined<-ScaleData(reyfmans.all.combined, verbose = FALSE)
reyfmans.all.combined<-RunPCA(reyfmans.all.combined, npcs = 30, verbose = FALSE)
reyfmans.all.combined<-RunUMAP(reyfmans.all.combined, reduction = "pca", dims = 1:18)
reyfmans.all.combined<-RunTSNE(reyfmans.all.combined, reduction = "pca", dims = 1:18)
reyfmans.all.combined<-FindNeighbors(reyfmans.all.combined, reduction = "pca", dims = 1:18)
reyfmans.all.combined<-FindClusters(reyfmans.all.combined, resolution = 0.1)
DimPlot(reyfmans.all.combined, reduction = "umap", label = TRUE)
exp_reyfmans<-AverageExpression(reyfmans.all.combined)
exp_reyfmans_rna<-exp_reyfmans[["RNA"]]
p1<-DimPlot(reyfmans.all.combined, reduction = "umap", group.by = "stim")
p2<-DimPlot(reyfmans.all.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(reyfmans.all.combined, reduction = "tsne", label = TRUE, repel = TRUE)
FeaturePlot(reyfmans.all.combined, "WDR31")
DefaultAssay(reyfmans.all.combined)<-"RNA"
```
**UMAP for TMEM231, IFT81 and NEK10**

![reyfman_gene_umap_1](https://user-images.githubusercontent.com/12661265/225610997-e5314613-d506-44d4-be1a-d8a4ef956902.png)


**Comparative Genomics (Phylogenetic profiling)**

**Load the following packages**

``` Python
library(data.table)
library(geneName)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(biomaRt)
library(ontologyIndex)
library(cluster)

newscores_50 <- filtering(pwd, 0.001, 50)
gower_dist_50 <- daisy(newscores_50[, 2:73], metric = "gower")
cluster_50 <- hclust(gower_dist_50, method = "ward.D2")
tree_50 <- cutree(cluster_50, k = 40)
aa_50 <- cbind(newscores_50[, c(1, 2)], tree_50) # Genes with clusters
aa_50_temp <- aa_50
colnames(aa_50_temp)[3] <- "cluster_num"
a <- as.data.frame(cbind(newscores_50, tree_50))
a[, 2:73] <- lapply(a[, 2:73], function(x) as.numeric(as.character(x)))
```



