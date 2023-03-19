**CilioGenics is an integrated method for predicting the ciliary genes**

To more accurately predict ciliary genes, CilioGenics combines several approaches, including as single-cell RNA sequencing, protein-protein interactions (PPIs), comparative genomics, transcription factor (TF)-network analysis, and text mining. 

The current repository contains codes for each method, and representative codes are presented below. 


**Single cell RNA-seq analysis**

**Load the following package**


``` Python
library(Seurat)
```

Dowload scRNA-seq file from NCBI and load the file.

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


**Annotate cells in scRNA-seq data**

``` Python
plot1 <- DimPlot(object = lungs.combined, reduction = "umap",label = TRUE, repel = FALSE,
        group.by = 'customclassif')
```        
![Umap_Cell_names](https://user-images.githubusercontent.com/12661265/225908354-d829eef2-8739-482b-8bd7-b044ecf21c16.png)

**UMAP for IFT88, TMEM231, NEK10, WDR54, WDR38, ZNF474, LGR5 (a non-ciliary gene)**


``` Python
FeaturePlot(lung.all.combined, features = c("IFT88", "TMEM231", "NEK10", "WDR54", 
                                            "WDR38", "ZNF474", "LGR5"),min.cutoff = "q10", 
            max.cutoff = "q90")
```    
![UMAP_IFT88](https://user-images.githubusercontent.com/12661265/225918548-1e3e476e-741d-467b-97bd-ca8562a402e2.png)



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



