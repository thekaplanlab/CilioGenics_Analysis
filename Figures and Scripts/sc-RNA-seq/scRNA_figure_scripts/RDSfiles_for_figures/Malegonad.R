
### theKaplanLab ##

### Hasan Can Demirci ##

## Required Libraries ##

library(Seurat)

library(biomaRt)

library(ggplot2)

## seting working directory ##

setwd("/home/kaplanlab/Desktop/HCD/sperm")

## importing raw rds file ##

sperm.combined <- readRDS("sperm.rds")

## finding variables ##

sperm.combined <- FindVariableFeatures(sperm.combined, selection.method = "vst", nfeatures = 2000)

## sacaling data ##

sperm.combined <- ScaleData(sperm.combined, verbose = FALSE)

## performing PCA to get 2D ##

sperm.combined <- RunPCA(sperm.combined, npcs = 30, verbose = FALSE)

## finding neighbors ##

sperm.combined <- FindNeighbors(sperm.combined, reduction = "pca", dims = 1:23)

## finding clusters ##

sperm.combined <- FindClusters(sperm.combined, resolution = 0.1) 

## UMAP arrangement ##

sperm.combined <- RunUMAP(sperm.combined, reduction = "pca", dims = 1:23)

## converting ensembel numbers to gene symbols ##

{
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  attributes <- c("ensembl_gene_id", "external_gene_name")
  ensembl_info <- getBM(attributes = attributes, 
                        filters = "ensembl_gene_id", 
                        values = rownames(sperm.combined), 
                        mart = ensembl)
  
  current_row_names <- rownames(GetAssayData(object = sperm.combined))
  
  new_row_names <- ensembl_info$external_gene_name[match(current_row_names, ensembl_info$ensembl_gene_id)]
  rownames(sperm.combined@assays$RNA@data) <- new_row_names
}

## naming Idents ##

sperm.combined <- RenameIdents(sperm.combined, 
                                `0` = "Sertoli Cells", 
                                `1` = "Mesenchymal cells", 
                                `2` = "Mesenchymal cells", 
                                `3` = "Mesenchymal cells", 
                                `4` = "Mesothelial cells", 
                                `5` = "Endothelial cells", 
                                `6` = "Germ cells", 
                                `7` = "Hematopoietic cells", 
                                `8` = "Epithelial cells", 
                                `9` = "Supporting cells", 
                                `10` = "Mural cells", 
                                `11` = "Leydig cells",
                                `12` = "Neural cells",
                                `13` ="Erythrocyte")

## getting umap plot ##

DimPlot(sperm.combined, reduction = "umap", label = TRUE)

## making large RDS as small RDS with DietSeurat function ##

sperm.combined <- DietSeurat(
  
  object = sperm.combined,
  counts = FALSE,
  data = TRUE,
  scale.data = FALSE,
  dimreducs = "umap"
  
)

## exporting new RDS to generate figures ##

saveRDS(sperm.combined, "sperm_dietseurat.RDS")


############################# Differential Analysis ############################

## finding markers ##

sperm_all_markers<-FindAllMarkers(sperm.combined, only.pos = TRUE, min.pct = 0.25,
                                  logfc.threshold = 0.25)


## converting ensembel IDs to gene symbols ##

{
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  attributes <- c("ensembl_gene_id", "external_gene_name")
  ensembl_info <- getBM(attributes = attributes, 
                        filters = "ensembl_gene_id", 
                        values = sperm_all_markers$gene
                        , 
                        mart = ensembl)
  
  current_row_names <- sperm_all_markers$gene
  
  
  new_row_names <- ensembl_info$external_gene_name[match(current_row_names, ensembl_info$ensembl_gene_id)]
  
  sperm_all_markers$gene_name <- new_row_names
}

## exporting markers table as CSV ##

write.csv(sperm_all_markers,"sperm_all_markers.csv")
