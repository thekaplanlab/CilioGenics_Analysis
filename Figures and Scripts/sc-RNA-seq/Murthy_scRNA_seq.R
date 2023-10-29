#Filtered data were downloaded
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178360 



setwd("/Users/oktayismailkaplan/Desktop/Manuscripts/Ciliogenics/Figures/Analysis")

library(Seurat)
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

#saveRDS(lungs.combined, "lungs.combined.RDS")

lungs.combined <- readRDS("lungs.combined.RDS", refhook = NULL)



lungs.combined<-ScaleData(lungs.combined, verbose = FALSE)
lungs.combined<-RunPCA(lungs.combined, npcs = 30, verbose = FALSE)
ElbowPlot(lungs.combined, ndims = 30)
lungs.combined<-RunUMAP(lungs.combined, reduction = "pca", dims = 1:23)
lungs.combined<-RunTSNE(lungs.combined, reduction = "pca", dims = 1:23)
lungs.combined<-FindNeighbors(lungs.combined, reduction = "pca", dims = 1:23)
lungs.combined<-FindClusters(lungs.combined, resolution = 0.1)

DimPlot(lungs.combined, reduction = "umap", label = TRUE)

########################################
#install.packages("HGNChelper") 
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)

library(dplyr)
library(HGNChelper)
library(openxlsx)
library(ggrepel)
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Lung" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = Murthy_PK[["integrated"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

cL_resutls = do.call("rbind", lapply(unique(Murthy_PK@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(Murthy_PK@meta.data[Murthy_PK@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(Murthy_PK@meta.data$seurat_clusters==cl)), 10)
}))

sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

Murthy_PK@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  Murthy_PK@meta.data$customclassif[Murthy_PK@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
library(ggrepel)
library(ggplot2)
plot1 <- DimPlot(object = Murthy_PK, reduction = "umap",label = TRUE, repel = FALSE,
        group.by = 'customclassif')
plot1 

library("dplyr")

FeaturePlot(object = lungs.combined, 
            features = c("IFT88", "IFT81", "WDR31", "WDR54", "ZNF474"),
            sort.cell = FALSE,
            min.cutoff = 'q2', 
            label = TRUE,
            repel = FALSE)


########################################

exp_lungs<-AverageExpression(lungs.combined)
exp_lungs_rna<-exp_lungs[["RNA"]]
exp_lungs_rna


#p1<-DimPlot(lungs.combined, reduction = "umap", group.by = "stim")
#p2<-DimPlot(lungs.combined, reduction = "umap", label = TRUE, repel = TRUE)
#p1 + p2

DimPlot(lungs.combined, reduction = "tsne", label = TRUE, repel = TRUE)
FeaturePlot(lungs.combined, "WDR31")

FeaturePlot(lungs.combined, "WDR54")

DefaultAssay(lungs.combined)<-"RNA"


##################################################

lungs.markers.rna <-FindAllMarkers(lungs.combined, only.pos = TRUE, min.pct = 0,
                                     logfc.threshold = 0)

lung_glist1<-lungs.markers.rna %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


lung.all.combined<-RenameIdents(lungs.combined, `0` = "Alveolar Macrophages", `1` = "Endothelial Cells", `2` = "Immune System Cells", 
                                    `3` = "Pulmonary alveolar type II cells", `4` = "Ciliated Cells", `5` = "Fibroblasts", `6` = "Basal cells (Airway progenitor cells)", `7` = "Fibroblasts", `8` = "Immune system cells", 
                                `9` = "Pulmonary alveolar type I cells", `10` = "Unknown", `11` = "Unknown", `12` = "Immune system cells",`13` = "Endothelial cell")


combined_lung_ciliated<-FindMarkers(lung.all.combined, ident.1 = "Ciliated Cells",
                                        min.pct = 0.1, logfc.threshold = 0)

View(combined_lung_ciliated)
combined_lung_ciliated$Gene_name<-rownames(combined_lung_ciliated)
combined_lung_ciliated$score<-combined_lung_ciliated$pct.1
combined_lung_ciliated$pct.2


ciliated <-lungs.markers.rna[lungs.markers.rna$cluster == 4,]
ciliated$gene

#search gene names
cilia_genes <- c("WDR54", "IFT88", "WDR31")

result <- filter(ciliated, grepl(paste(cilia_genes, collapse="|"), gene))
result

#grep(pattern = 'WDR54|IFT88', x = ciliated$gene)

#match(c("WDR54", "IFT88", "WDR31"), ciliated$gene)

#lungs_reduced <-DietSeurat(lung.all.combined, counts = TRUE, 
                           #assays = "integrated", dimreducs = "umap")

lungs_reduced1 <-DietSeurat(lung.all.combined, counts = TRUE, 
                           assays = "RNA", dimreducs = "umap")

DimPlot(lungs_reduced1, reduction = "umap", label = TRUE, repel = FALSE)

VlnPlot(lung.all.combined, features = c("WDR54", "WDR31"), slot = "counts", log = TRUE)

#FeaturePlot(lung.all.combined, features = c("IFT88", "NEK10", "WDR38", "IFT22", 
                                           # "TMEM231", "IFT27", "ZNF474", "ZC2HC1A", 
                                            #"ARMC9"))

FeaturePlot(lung.all.combined, features = c("IFT88", "TMEM231", "NEK10", "WDR54", 
                                            "WDR38", "ZNF474", "LGR5"),min.cutoff = "q10", 
            max.cutoff = "q90")

FeaturePlot(lung.all.combined, features = c("IFT88", "TMEM231", "NEK10", "WDR54", 
                                            "WDR38", "LGR5"),min.cutoff = "q10", 
            max.cutoff = "q90")

#DotPlot(lung.all.combined, features = c("IFT88", "TMEM231", "NEK10", "WDR54", 
                                        #"WDR38", "LGR5"))

library(ggplot2)
DotPlot(lung.all.combined, features = c("IFT88", "TMEM231", "NEK10", "WDR54", 
                                        "WDR38", "LGR5"),
        cols = c("blue","red"), 
        dot.scale = 4)+
  RotatedAxis()+
  coord_flip()+
  theme(axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7)) + 
  xlab('Cluster') +  
  ylab('Gene')
  

DotPlot(lung.all.combined, features = c("NEK10", "WDR54", 
                                        "WDR38", "LGR5","C1orf194", "TMEM190", "WDR31", "ZC2HC1A", "IFT88", "TMEM231"),
        cols = c("blue","red"), 
        dot.scale = 4)+
  RotatedAxis()+
  coord_flip()+
  theme(axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7)) + 
  xlab('Cluster') +  
  ylab('Gene')





  #Single cell heatmap of feature expression
DoHeatmap(subset(lung.all.combined, downsample = 100), features = features, size = 3)

#saveRDS(lung.all.combined, "lung_seurat.RDS")
#lung.all.combined<-readRDS("reyfmans_seurat.RDS")
#saveRDS(lungs_reduced, "lung_seurat_reduced.RDS")
