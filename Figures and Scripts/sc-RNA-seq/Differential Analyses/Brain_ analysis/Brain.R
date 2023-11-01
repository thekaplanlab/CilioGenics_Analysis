
### theKaplanLab ##

### Hasan Can Demirci ##

## Required Libraries ##

library(Seurat)

library(biomaRt)

library(ggplot2)

## seting working directory ##

setwd("/home/kaplanlab/Desktop/HCD/Brain")

## importing raw rds file ##

brain.combined <- readRDS("brain.rds")

## finding variables ##

brain.combined <- FindVariableFeatures(brain.combined, selection.method = "vst", nfeatures = 2000)

## sacaling data ##

brain.combined <- ScaleData(brain.combined, verbose = FALSE)

## performing PCA to get 2D ##

brain.combined <- RunPCA(brain.combined, npcs = 30, verbose = FALSE)

## finding neighbors ##

brain.combined <- FindNeighbors(brain.combined, reduction = "pca", dims = 1:23)

## finding clusters ##

brain.combined <- FindClusters(brain.combined, resolution = 0.1) 

## UMAP arrangement ##

brain.combined <- RunUMAP(brain.combined, reduction = "pca", dims = 1:23)

## converting ensembel numbers to gene symbols ##

{
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  attributes <- c("ensembl_gene_id", "external_gene_name")
  ensembl_info <- getBM(attributes = attributes, 
                        filters = "ensembl_gene_id", 
                        values = rownames(brain.combined), 
                        mart = ensembl)
  
  current_row_names <- rownames(GetAssayData(object = brain.combined))
  
  new_row_names <- ensembl_info$external_gene_name[match(current_row_names, ensembl_info$ensembl_gene_id)]
  rownames(brain.combined@assays$RNA@data) <- new_row_names
}

## naming Idents ##

brain.combined <- RenameIdents(brain.combined, 
                               `0` = "Neurons", 
                               `1` = "Neurons", 
                               `2` = "Oligodendrocytes", 
                               `3` = "Neurons", 
                               `4` = "Bergmann glia cells", 
                               `5` = "Oligodendrocytes", 
                               `6` = "Immune cells", 
                               `7` = "Neurons", 
                               `8` = "Oligodendrocyte precursor cells", 
                               `9` = "Vascular cells", 
                               `10` = "Neurons") 



## getting umap plot ##

DimPlot(brain.combined, reduction = "umap",label = TRUE)

## making large RDS as small RDS with DietSeurat function ##

brain.combined <- DietSeurat(
  
  object = brain.combined,
  counts = FALSE,
  data = TRUE,
  scale.data = FALSE,
  dimreducs = "umap"
  
)

## exporting new RDS to generate figures ##

saveRDS(brain.combined, "brain_dietseurat.RDS")

############################# Differential Analysis ############################

## finding markers ##

brain_all_markers<-FindAllMarkers(brain.combined, only.pos = TRUE, min.pct = 0.25,
                                  logfc.threshold = 0.25)


## converting ensembel IDs to gene symbols ##

{
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  attributes <- c("ensembl_gene_id", "external_gene_name")
  ensembl_info <- getBM(attributes = attributes, 
                        filters = "ensembl_gene_id", 
                        values = brain_all_markers$gene
                        , 
                        mart = ensembl)
  
  current_row_names <- brain_all_markers$gene
  
  
  new_row_names <- ensembl_info$external_gene_name[match(current_row_names, ensembl_info$ensembl_gene_id)]
  
  brain_all_markers$gene_name <- new_row_names
}

## exporting markers table as CSV ##

write.csv(brain_all_markers,"brain_all_markers.csv")





