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
lung_dir <- dir("./GSE178360_RAW")
lung_dir_1 <- lung_dir[grepl("filtered", lung_dir)]
lungs <- list()
#install.packages("hdf5r") if it is not installed.
for (i in 1:3){
  lung_mat <- Read10X_h5(lung_dir_1[i])
  lung_mat1 <- CreateSeuratObject(lung_mat)
  lungs[[i]]<-lung_mat1
}
lungs.list<-lapply(X = lungs, FUN = function(x) {
  x<-NormalizeData(x)
  x<-FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
#ctr+shft+alt+m replace all
features<-SelectIntegrationFeatures(object.list = lungs.list)
lungs.anchors<-FindIntegrationAnchors(object.list = lungs.list, anchor.features = features)
integrate_all <- Reduce(intersect, lapply(lungs.anchors@object.list, rownames)) # create list of common genes to keep
lungs.combined<-IntegrateData(anchorset = lungs.anchors, features.to.integrate = integrate_all)
DefaultAssay(lungs.combined)<-"integrated"
```

**Data visulaziation**
``` Python
lungs.combined <- readRDS("lungs.combined.RDS", refhook = NULL)
lungs.combined<-ScaleData(lungs.combined, verbose = FALSE)
lungs.combined<-RunPCA(lungs.combined, npcs = 30, verbose = FALSE)
ElbowPlot(lungs.combined, ndims = 30)
lungs.combined<-RunUMAP(lungs.combined, reduction = "pca", dims = 1:23)
lungs.combined<-RunTSNE(lungs.combined, reduction = "pca", dims = 1:23)
lungs.combined<-FindNeighbors(lungs.combined, reduction = "pca", dims = 1:23)
lungs.combined<-FindClusters(lungs.combined, resolution = 0.1)
DimPlot(lungs.combined, reduction = "umap", label = TRUE)
```
**UMAP for TMEM231, IFT81 and NEK10**

![reyfman_gene_umap_1](https://user-images.githubusercontent.com/12661265/225610997-e5314613-d506-44d4-be1a-d8a4ef956902.png)

**Annotate cells in scRNA-seq data**

``` Python
plot1 <- DimPlot(object = lungs.combined, reduction = "umap",label = TRUE, repel = FALSE,
        group.by = 'customclassif')
```        

![UMAP_names](https://user-images.githubusercontent.com/12661265/225907931-ee5e83dc-964f-46c5-a277-09e7a21f3779.png)

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



