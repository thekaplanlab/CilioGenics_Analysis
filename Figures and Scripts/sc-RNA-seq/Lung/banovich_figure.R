library(Seurat)
library(ggplot2)

a <- readRDS("/home/kaplanlab/Desktop/Desktop_Items/MUSTAFA_PIR/web/CilioGenics/data/banovich_reduced.RDS")

DotPlot(a, features = c("WDR38", "WDR54", "ZC2HC1A", "TMEM190", "WDR31", "ZNF474"))

ggsave("banovich-dotplot-new.tiff", units="in", bg="white" ,width= 11, height = 7, dpi = 300)


DimPlot(a, reduction = "umap", label=TRUE)

ggsave("banovich-dimplot-new.tiff", units="in", width = 9  ,height = 7, dpi = 300)


FeaturePlot(object = a, 
            features = c("WDR38", "WDR54", "ZC2HC1A", "TMEM190", "WDR31", "ZNF474"))

ggsave("banovich-candidate-new.tiff", units="in", width = 9  ,height = 7, dpi = 300)


FeaturePlot(object = a, 
            features = c("IFT88", "TMEM231", "NEK10"))

ggsave("banovich-markers-new.tiff", units="in", width = 9  ,height = 7, dpi = 300)
