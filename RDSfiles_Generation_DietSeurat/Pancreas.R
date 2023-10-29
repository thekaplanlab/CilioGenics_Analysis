
### theKaplanLab ##

### Hasan Can Demirci ##

## Required Libraries ##

library(Seurat)

library(biomaRt)

library(ggplot2)

## seting working directory ##

setwd("/home/kaplanlab/Desktop/HCD/pancreas")

## importing raw rds file ##

pancreas.combined <- readRDS("pancreas.rds")

## finding variables ##

pancreas.combined <- FindVariableFeatures(pancreas.combined, selection.method = "vst", nfeatures = 2000)

## sacaling data ##

pancreas.combined <- ScaleData(pancreas.combined, verbose = FALSE)

## performing PCA to get 2D ##

pancreas.combined <- RunPCA(pancreas.combined, npcs = 30, verbose = FALSE)

## finding neighbors ##

pancreas.combined <- FindNeighbors(pancreas.combined, reduction = "pca", dims = 1:23)

## finding clusters ##

pancreas.combined <- FindClusters(pancreas.combined, resolution = 0.1) 

## UMAP arrangement ##

pancreas.combined <- RunUMAP(pancreas.combined, reduction = "pca", dims = 1:23)

## converting ensembel numbers to gene symbols ##

{
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  attributes <- c("ensembl_gene_id", "external_gene_name")
  ensembl_info <- getBM(attributes = attributes, 
                        filters = "ensembl_gene_id", 
                        values = rownames(pancreas.combined), 
                        mart = ensembl)
  
  current_row_names <- rownames(GetAssayData(object = pancreas.combined))
  
  new_row_names <- ensembl_info$external_gene_name[match(current_row_names, ensembl_info$ensembl_gene_id)]
  rownames(pancreas.combined@assays$RNA@data) <- new_row_names
}

## naming Idents ##

pancreas.combined <- RenameIdents(pancreas.combined, 
                                `0` = "Type A enteroendocrine cells", 
                                `1` = "Pancreatic ductal cells", 
                                `2` = "Acinar cell", 
                                `3` = "Type B pancreatic cell types", 
                                `4` = "Type A enteroendocrine cells", 
                                `5` = "Native cells", 
                                `6` = "Mesenchymal cells", 
                                `7` = "type D enteroendocrine cell") 
                              
## getting umap plot ##

DimPlot(pancreas.combined, reduction = "umap",label = TRUE)

## making large RDS as small RDS with DietSeurat function ##

pancreas.combined <- DietSeurat(
  
  object = pancreas.combined,
  counts = FALSE,
  data = TRUE,
  scale.data = FALSE,
  dimreducs = "umap"
  
)

## exporting new RDS to generate figures ##

saveRDS(pancreas.combined, "pancreas_dietseurat.RDS")

############################# Differential Analysis ############################

## finding markers ##

pancreas_all_markers<-FindAllMarkers(pancreas.combined, only.pos = TRUE, min.pct = 0.25,
                                  logfc.threshold = 0.25)


## converting ensembel IDs to gene symbols ##

{
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  attributes <- c("ensembl_gene_id", "external_gene_name")
  ensembl_info <- getBM(attributes = attributes, 
                        filters = "ensembl_gene_id", 
                        values = pancreas_all_markers$gene
                        , 
                        mart = ensembl)
  
  current_row_names <- pancreas_all_markers$gene
  
  
  new_row_names <- ensembl_info$external_gene_name[match(current_row_names, ensembl_info$ensembl_gene_id)]
  
  pancreas_all_markers$gene_name <- new_row_names
}

## exporting markers table as CSV ##

write.csv(pancreas_all_markers,"pancreas_all_markers.csv")