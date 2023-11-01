library(Seurat)
library(ggplot2)

a <- readRDS("~/Desktop/HCD/revision/liver_dietseurat/liver_dietseurat.rds")

DotPlot(a, features = c("WDR38", "WDR54", "ZC2HC1A", "TMEM190", "WDR31", "ZNF474"))

ggsave("liver-dotplot-new.tiff", units="in", bg="white" ,width= 11, height = 7, dpi = 300)


DimPlot(a, reduction = "umap", label=TRUE)

ggsave("liver-dimplot-new.tiff", units="in", width = 9  ,height = 7, dpi = 300)


FeaturePlot(object = a, 
            features = c("WDR38", "WDR54", "ZC2HC1A", "TMEM190", "WDR31", "ZNF474"))

ggsave("liver-candidate-new.tiff", units="in", width = 9  ,height = 7, dpi = 300)


FeaturePlot(object = a, 
            features = c("IFT88", "TMEM231", "NEK10"))

ggsave("liver-markers-new.tiff", units="in", width = 9  ,height = 7, dpi = 300)
