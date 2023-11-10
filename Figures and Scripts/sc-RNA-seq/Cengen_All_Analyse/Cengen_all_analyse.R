###The KaplanLab###
###Advisor: Oktay I. Kaplan###
###Author: Mustafa Emre KORKMAZ###

###Publication: Hammarlund, M., Hobert, O., Miller, D. M., & Sestan, N. (2018). The CeNGEN project: the complete gene expression map of an entire nervous system --- https://doi.org/10.1016/j.neuron.2018.07.042 
###rds_file: http://storage.cengen.org/100720_L4_neuron_only_Seurat_with_RIV_SMD_stressed.rds -- C. elegans Neurons (70,296 cells)

###Necessary libraries are loaded.
library(Seurat)
library(biomaRt)
library(ggplot2)
library(dplyr)

###Reading a dataset from an RDS file.
celegans.neuron <- readRDS("100720_L4_neuron_only_Seurat_with_RIV_SMD_stressed.rds")

###Selecting important features (genes) from the data.
celegans.neuron <- FindVariableFeatures(celegans.neuron, selection.method = "vst", nfeatures = 2000)

###Scaling the data for analysis.
celegans.neuron <- ScaleData(celegans.neuron, verbose = FALSE)

###Reducing dimensionality using PCA.
celegans.neuron <- RunPCA(celegans.neuron, npcs = 30, verbose = FALSE)

###Calculating data point neighborhoods.
celegans.neuron <- FindNeighbors(celegans.neuron, reduction = "pca", dims = 1:23)

###Clustering data points based on neighborhoods.
celegans.neuron <- FindClusters(celegans.neuron, resolution = 0.1) 

###Reducing dimensionality using UMAP (Uniform Manifold Approximation and Projection).
celegans.neuron <- RunUMAP(celegans.neuron, reduction = "pca", dims = 1:23)

###Assigning identities or labels to clusters.
Idents(celegans.neuron) <- celegans.neuron@meta.data$Cell.type

###Creating a UMAP-based plot for data visualization.
DimPlot(celegans.neuron, reduction = "umap",label = TRUE)

###Read csv file to get homology list
homology_list <- read.csv("C_elegans_Human_Homologs_20.01.2019 - UPDATED.csv")

###Converting Wormbase IDs into Human symbol in RDS data.
new_row_names <- ifelse(rownames(celegans.neuron) %in% homology_list$WBGeneID,
                        homology_list$Human_Symbol[match(rownames(celegans.neuron),homology_list$WBGeneID)],
                        "")
rownames(celegans.neuron@assays$SCT@data) <- new_row_names

###Identifying markers (all genes) for each cluster. 
celegans.neuron.markers.neww <- FindAllMarkers(celegans.neuron ,only.pos = TRUE, logfc.threshold = 0.01, min.pct = 0.01,test.use = "wilcox")

###Give a score for expression of each gene
rownames(celegans.neuron)<-rownames(celegans.neuron.markers.neww)
celegans.neuron.markers.neww$score<-celegans.neuron.markers.neww$pct.1-
  celegans.neuron.markers.neww$pct.2

###Converting Wormbase IDs into Human symbol.
celegans.neuron.markers.neww$human_symbol <- ifelse(celegans.neuron.markers.neww$gene %in% homology_list$WBGeneID,
                                                    homology_list$Human_Symbol[match(celegans.neuron.markers.neww$gene,homology_list$WBGeneID)],
                                                    "")

write.csv(celegans.neuron.markers.neww,"cengen_all_scores_human.csv")

###Easy Threshold Gene List Check if Idents(celegans.neuron) <- celegans.neuron@meta.data$modality
celegans.sensory.0025 <- celegans.neuron.markers.neww %>%
  filter(pct.1>=0.025 & cluster=="Sensory")

#Feature plot for marker genes
FeaturePlot(object = celegans.neuron, 
            features = c("IFT88", "TMEM231","NEK10"),
            sort.cell = TRUE,
            min.cutoff = 'q2', 
            label = TRUE,
            repel = TRUE)

#Feature plot for candidate genes
FeaturePlot(object = celegans.neuron, 
            features = c("WDR38", "TMEM145","C1orf194","TMEM190","ZC2HC1A","WDR54", "ZNF474","WDR31"),
            sort.cell = TRUE,
            min.cutoff = NA, 
            label = TRUE,
            repel = TRUE)

#Dot plot for candidate genes
DotPlot(a, features = c("WDR38", "WDR54", "ZC2HC1A", "TMEM190", "WDR31", "ZNF474"))

########################################################
cengen_diet <- DietSeurat(celegans.neuron,
                          counts = TRUE,
                          data = TRUE,
                          scale.data = FALSE,
                          features = NULL,
                          assays = NULL,
                          dimreducs = NULL,
                          graphs = NULL,
                          misc = TRUE)

saveRDS(cengen_diet, "cengen_dietseurat.rds")
############################################################