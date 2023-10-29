
### theKaplanLab ##

### Hasan Can Demirci ##

## Required Libraries ##

library(Seurat)

library(biomaRt)

library(ggplot2)

## seting working directory ##

setwd("/home/kaplanlab/Desktop/HCD/trachea")

## importing raw rds file ##

trachea.combined <- readRDS("trachea.rds")

## finding variables ##

trachea.combined <- FindVariableFeatures(trachea.combined, selection.method = "vst", nfeatures = 2000)

## sacaling data ##

trachea.combined <- ScaleData(trachea.combined, verbose = FALSE)

## performing PCA to get 2D ##

trachea.combined <- RunPCA(trachea.combined, npcs = 30, verbose = FALSE)

## finding neighbors ##

trachea.combined <- FindNeighbors(trachea.combined, reduction = "pca", dims = 1:23)

## finding clusters ##

trachea.combined <- FindClusters(trachea.combined, resolution = 0.1) 

## UMAP arrangement ##

trachea.combined <- RunUMAP(trachea.combined, reduction = "pca", dims = 1:23)

## converting ensembel numbers to gene symbols ##

{
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  attributes <- c("ensembl_gene_id", "external_gene_name")
  ensembl_info <- getBM(attributes = attributes, 
                        filters = "ensembl_gene_id", 
                        values = rownames(trachea.combined), 
                        mart = ensembl)
  
  current_row_names <- rownames(GetAssayData(object = trachea.combined))
  
  new_row_names <- ensembl_info$external_gene_name[match(current_row_names, ensembl_info$ensembl_gene_id)]
  rownames(trachea.combined@assays$RNA@data) <- new_row_names
}

## naming Idents ##

trachea.combined <- RenameIdents(trachea.combined, 
                                `0` = "Basal Cells", 
                                `1` = "Connective Tissue Cells \n Fibroblasts", 
                                `2` = "Neutrophils", 
                                `3` = "CD8-positive T Cells", 
                                `4` = "Macrophages", 
                                `5` = "B Cells", 
                                `6` = "Goblet Cells", 
                                `7` = "Ciliated Cells",
                                `8` = "Endothelial Cells",
                                `9` = "Basal Cells",
                                `10`= "Plasma Cells",
                                `11` = "Smooth Muscles Cells",
                                `12` = "Basal Cells",
                                `13` = "Mast Cells") 


## getting umap plot ##

DimPlot(trachea.combined, reduction = "umap",label = TRUE)

## making large RDS as small RDS with DietSeurat function ##

trachea.combined <- DietSeurat(
  
  object = trachea.combined,
  counts = FALSE,
  data = TRUE,
  scale.data = FALSE,
  dimreducs = "umap"
  
)

## exporting new RDS to generate figures ##

saveRDS(trachea.combined, "trachea_dietseurat.RDS")

############################# Differential Analysis ############################

## finding markers ##

trachea_all_markers<-FindAllMarkers(trachea.combined, only.pos = TRUE, min.pct = 0.25,
                                  logfc.threshold = 0.25)


## converting ensembel IDs to gene symbols ##

{
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  attributes <- c("ensembl_gene_id", "external_gene_name")
  ensembl_info <- getBM(attributes = attributes, 
                        filters = "ensembl_gene_id", 
                        values = trachea_all_markers$gene
                        , 
                        mart = ensembl)
  
  current_row_names <- trachea_all_markers$gene
  
  
  new_row_names <- ensembl_info$external_gene_name[match(current_row_names, ensembl_info$ensembl_gene_id)]
  
  trachea_all_markers$gene_name <- new_row_names
}

## exporting markers table as CSV ##

write.csv(trachea_all_markers,"trachea_all_markers.csv")








