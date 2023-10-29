
### theKaplanLab ##

### Hasan Can Demirci ##

## Required Libraries ##

library(Seurat)

library(biomaRt)

library(ggplot2)

## seting working directory ##

setwd("/home/kaplanlab/Desktop/HCD/retina")

## importing raw rds file ##

retina.combined <- readRDS("retina.rds")

## finding variables ##

retina.combined <- FindVariableFeatures(retina.combined, selection.method = "vst", nfeatures = 2000)

## sacaling data ##

retina.combined <- ScaleData(retina.combined, verbose = FALSE)

## performing PCA to get 2D ##

retina.combined <- RunPCA(retina.combined, npcs = 30, verbose = FALSE)

## finding neighbors ##

retina.combined <- FindNeighbors(retina.combined, reduction = "pca", dims = 1:23)

## finding clusters ##

retina.combined <- FindClusters(retina.combined, resolution = 0.1) 

## UMAP arrangement ##

retina.combined <- RunUMAP(retina.combined, reduction = "pca", dims = 1:23)

## converting ensembel numbers to gene symbols ##

{
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  attributes <- c("ensembl_gene_id", "external_gene_name")
  ensembl_info <- getBM(attributes = attributes, 
                        filters = "ensembl_gene_id", 
                        values = rownames(retina.combined), 
                        mart = ensembl)
  
  current_row_names <- rownames(GetAssayData(object = retina.combined))
  
  new_row_names <- ensembl_info$external_gene_name[match(current_row_names, ensembl_info$ensembl_gene_id)]
  rownames(retina.combined@assays$RNA@data) <- new_row_names
}

## naming Idents ##

retina.combined <- RenameIdents(retina.combined, 
                                `0` = "Retinal rod cells", 
                                `1` = "Retinal rod cells", 
                                `2` = "Retinal rod cells", 
                                `3` = "Mueller cells", 
                                `4` = "OFF−bipolar cells", 
                                `5` = "Retinal cone cells", 
                                `6` = "ON−bipolar cells", 
                                `7` = "Retinal ganglion cells", 
                                `8` = "ON−bipolar cells", 
                                `9` = "Microglial cells", 
                                `10` = "Amacrine cells", 
                                `11` = "Native cells")

## getting umap plot ##

 DimPlot(retina.combined, reduction = "umap", label = TRUE)
 
 ## making large RDS as small RDS with DietSeurat function ##
 
 retina.combined <- DietSeurat(
   
   object = retina.combined,
   counts = FALSE,
   data = TRUE,
   scale.data = FALSE,
   dimreducs = "umap"
   
 )
 
 ## exporting new RDS to generate figures ##
 
 saveRDS(retina.combined, "retina_dietseurat.RDS")
 
 ############################# Differential Analysis ############################
 
 ## finding markers ##
 
 retina_all_markers<-FindAllMarkers(retina.combined, only.pos = TRUE, min.pct = 0.25,
                                   logfc.threshold = 0.25)
 
 
 ## converting ensembel IDs to gene symbols ##
 
 {
   ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
   attributes <- c("ensembl_gene_id", "external_gene_name")
   ensembl_info <- getBM(attributes = attributes, 
                         filters = "ensembl_gene_id", 
                         values = retina_all_markers$gene
                         , 
                         mart = ensembl)
   
   current_row_names <- retina_all_markers$gene
   
   
   new_row_names <- ensembl_info$external_gene_name[match(current_row_names, ensembl_info$ensembl_gene_id)]
   
   retina_all_markers$gene_name <- new_row_names
 }
 
 ## exporting markers table as CSV ##
 
 write.csv(retina_all_markers,"retina_all_markers.csv")
 
 
 
 
 
 

