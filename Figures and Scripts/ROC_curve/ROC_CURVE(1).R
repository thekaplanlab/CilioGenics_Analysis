
## Required libraries ##

library(data.table)
library(pROC)
library(purrr)
library(dplyr)
library(ggplot2)

## working direcory ##

setwd("C:\\Users\\Hasan Can Demirci\\Desktop\\CilioGenics Revision\\ROC_Curve")

## data import ##

all_scores1<-fread("score_table_roc.txt")

## roc creation ## 

rocs <- roc(is.ciliary ~ interaction_scores + phylogenetic_scores + motif_score +
              pa_score + single_cell_score + ciliacarta + total_score_sums, data = all_scores1)

## labels naming ##

lbls<-c("Interactions", "Comparative genomics", "TF binding", "Text mining", "scRNA-seq", "CiliaCarta", "CilioGenics")

## AUC calculation ##
rocs %>% 
  map(~tibble(AUC = .x$auc)) %>% 
  rbindlist() -> data.auc

data.auc$Methods<-lbls

## new label creating ##
data.auc %>% 
  mutate(label_long=paste0(Methods," , AUC = ",paste(round(AUC,2))),
         label_AUC=paste0("AUC = ",paste(round(AUC,2)))) -> data.labels

## get plot ##

ggroc(rocs, size = 4) +
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#CC79A7")) +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1),
               color="darkgrey", linetype="dashed") +
  labs(color='Methods') + theme_bw(base_size = 30) +
  theme(legend.key.size = unit(2, "cm"))

## export plot ##

ggsave("roc_all_bw_ccarta-xx-3.jpg", width = 21, height = 14, dpi = 300)

