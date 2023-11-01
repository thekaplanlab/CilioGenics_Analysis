
### theKaplanLab ##

### Hasan Can Demirci ##

## Required Libraries ##

library(Seurat)

library(biomaRt)

library(ggplot2)

## seting working directory ##

setwd("/home/kaplanlab/Desktop/HCD/liver")

## importing raw rds file ##

liver.combined <- readRDS("liver.rds")

## finding variables ##

liver.combined <- FindVariableFeatures(liver.combined, selection.method = "vst", nfeatures = 2000)

## sacaling data ##

liver.combined <- ScaleData(liver.combined, verbose = FALSE)

## performing PCA to get 2D ##

liver.combined <- RunPCA(liver.combined, npcs = 30, verbose = FALSE)

## finding neighbors ##

liver.combined <- FindNeighbors(liver.combined, reduction = "pca", dims = 1:23)

## finding clusters ##

liver.combined <- FindClusters(liver.combined, resolution = 0.1) 

## UMAP arrangement ##

liver.combined <- RunUMAP(liver.combined, reduction = "pca", dims = 1:23)

## converting ensembel numbers to gene symbols ##

{
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  attributes <- c("ensembl_gene_id", "external_gene_name")
  ensembl_info <- getBM(attributes = attributes, 
                        filters = "ensembl_gene_id", 
                        values = rownames(liver.combined), 
                        mart = ensembl)
  
  current_row_names <- rownames(GetAssayData(object = liver.combined))
  
  new_row_names <- ensembl_info$external_gene_name[match(current_row_names, ensembl_info$ensembl_gene_id)]
  rownames(liver.combined@assays$RNA@data) <- new_row_names
}

## naming Idents ##

liver.combined <- RenameIdents(liver.combined, 
                               `0` = "Immune cells", 
                               `1` = "Hepatocytes", 
                               `2` = "Hepatocytes", 
                               `3` = "Endotehelial cell of hepatic sinusoids", 
                               `4` = "Immune cells", 
                               `5` = "Endothelial cells", 
                               `6` = "Neutrophils", 
                               `7` = "Erythrocytes", 
                               `8` = "Fibroblasts", 
                               `9` = "Intrahepatic Cholangiocyte", 
                               `10` = "Immune cells") 
                               
## getting umap plot ##

DimPlot(liver.combined, reduction = "umap", label = TRUE)

## making large RDS as small RDS with DietSeurat function ##

liver.combined <- DietSeurat(
  
  object = liver.combined,
  counts = FALSE,
  data = TRUE,
  scale.data = FALSE,
  dimreducs = "umap"
  
)

## exporting new RDS to generate figures ##

saveRDS(liver.combined, "liver_dietseurat.RDS")


############################# Differential Analysis ############################

## finding markers ##

liver_all_markers<-FindAllMarkers(liver.combined, only.pos = TRUE, min.pct = 0.25,
                                  logfc.threshold = 0.25)


## converting ensembel IDs to gene symbols ##

{
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  attributes <- c("ensembl_gene_id", "external_gene_name")
  ensembl_info <- getBM(attributes = attributes, 
                        filters = "ensembl_gene_id", 
                        values = liver_all_markers$gene
                        , 
                        mart = ensembl)
  
  current_row_names <- liver_all_markers$gene
  
  
  new_row_names <- ensembl_info$external_gene_name[match(current_row_names, ensembl_info$ensembl_gene_id)]
  
  liver_all_markers$gene_name <- new_row_names
}

## exporting markers table as CSV ##

write.csv(liver_all_markers,"liver_all_markers.csv")

