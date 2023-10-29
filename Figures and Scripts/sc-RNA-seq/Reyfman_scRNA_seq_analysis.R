

library(Seurat)

reyfman_dir<-dir("./sc_files/Reyfman")

reyfman_dir<-reyfman_dir[grepl("filtered", reyfman_dir)]
reyfman_dir_don<-reyfman_dir[grepl("Donor", reyfman_dir)]
reyfman_dir_dis<-reyfman_dir[!grepl("Donor", reyfman_dir)]
reyfman_dir_dis<-reyfman_dir_dis[!grepl("Cryobiopsy", reyfman_dir_dis)]

reyfman_all<-c(reyfman_dir_don, reyfman_dir_dis)

reyfmans<-list()

for (i in 1:8){
  reyfman<-Read10X_h5(paste0("sc_files/Reyfman/", reyfman_all[i]))
  reyfmans1<-CreateSeuratObject(reyfman)
  reyfmans1$samples<-paste0("Donor_",i)
  reyfmans1$stim<-"control"
  reyfmans[[i]]<-reyfmans1
}

reyfmans_d<-list()
for (i in 9:16){
  reyfman<-Read10X_h5(paste0("sc_files/Reyfman/", reyfman_all[i]))
  reyfmans1<-CreateSeuratObject(reyfman)
  reyfmans1$samples<-paste0("Dis_",i)
  reyfmans1$stim<-"disease"
  reyfmans_d[[i-8]]<-reyfmans1
}

reyfmans1.list<-lapply(X = reyfmans, FUN = function(x) {
  x<-NormalizeData(x)
  x<-FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

reyfmans2.list<-lapply(X = reyfmans_d, FUN = function(x) {
  x<-NormalizeData(x)
  x<-FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features<-SelectIntegrationFeatures(object.list = reyfmans1.list)
reyfmans1.anchors<-FindIntegrationAnchors(object.list = reyfmans1.list, anchor.features = features)
reyfmans1.combined<-IntegrateData(anchorset = reyfmans1.anchors)

features2<-SelectIntegrationFeatures(object.list = reyfmans2.list)
reyfmans2.anchors<-FindIntegrationAnchors(object.list = reyfmans2.list, anchor.features = features2)
reyfmans2.combined<-IntegrateData(anchorset = reyfmans2.anchors)

reyfmans_all_list<-list(reyfmans1.combined, reyfmans2.combined)

features3<-SelectIntegrationFeatures(object.list = reyfmans_all_list)
reyfmans.all.anchors<-FindIntegrationAnchors(object.list = reyfmans_all_list, anchor.features = features3)
reyfmans.all.combined<-IntegrateData(anchorset = reyfmans.all.anchors)

DefaultAssay(reyfmans.all.combined)<-"integrated"

#saveRDS(reyfmans.all.combined, "reyfmans.all.combined2.RDS")

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

# find markers for every cluster compared to all remaining cells, report only the positive ones
reyfmans.markers.rna<-FindAllMarkers(reyfmans.all.combined, only.pos = TRUE, min.pct = 0,
                                     logfc.threshold = 0)

glist1<-reyfmans.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


reyfmans.all.combined<-RenameIdents(reyfmans.all.combined, `0` = "Macrophages", `1` = "AT2 Cells", `2` = "Kupffer Cells", 
                                    `3` = "AT1 Cells", `4` = "Ciliated Cells", `5` = "T Cells", `6` = "B Cells", `7` = "7", `8` = "Endothelial Cells", `9` = "9", 
                                    `10` = "Granulocytes", `11` = "Fibroblasts", `12` = "12")

combined_reyfmans_ciliated<-FindMarkers(reyfmans.all.combined, ident.1 = "Ciliated Cells",
                                        min.pct = 0.1, logfc.threshold = 0)
combined_reyfmans_ciliated$Gene_name<-rownames(combined_reyfmans_ciliated)
combined_reyfmans_ciliated$score<-combined_reyfmans_ciliated_min$pct.1-
  combined_reyfmans_ciliated_min$pct.2


ciliated<-reyfmans.markers[reyfmans.markers$cluster == 4,]
reyfmans.reduced<-DietSeurat(reyfmans.all.combined, counts = FALSE, assays = "RNA", dimreducs = "umap")
DimPlot(reyfmans.reduced, reduction = "umap", label = TRUE, repel = TRUE)


# saveRDS(reyfmans.all.combined, "reyfmans_seurat.RDS")
# reyfmans.all.combined<-readRDS("reyfmans_seurat.RDS")
# saveRDS(reyfmans.reduced, "reyfmans_seurat_reduced.RDS")

# reyfman_data<-readRDS("files/reyfmans_seurat.RDS")
# reyfman<-readRDS("files/markers_reyfman.RDS")
