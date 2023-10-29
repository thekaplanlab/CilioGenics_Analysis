
source("generate_objects.R")

## Single cell scores ----

### Reyfman ----

reyfman <- readRDS("files/markers_reyfman.RDS") %>%
  filter(
    cluster == "Ciliated Cells", pct.1 - pct.2 > 0,
    p_val_adj < 0.001
  ) %>%
  hgncConverter("gene") %>%
  mutate(
    mean_score = normalization0(log(pct.1 - pct.2)),
    source = "reyfman"
  ) %>%
  rename("Gene_name" = "gene") %>%
  dplyr::select(Gene_name, mean_score, source) %>%
  group_by(Gene_name) %>%
  mutate(mean_score = max(mean_score)) %>%
  distinct()
reyfman$mean_score[reyfman$mean_score == 0] <- 0.00001


### Carraro et al. ----

carraro <- readRDS("files/markers.RDS")
carraro <- carraro[[2]] %>%
  filter(p_val_adj < 0.001) %>%
  mutate(
    mean_score = pct.1 - pct.2,
    source = "carraro"
  ) %>%
  tibble::rownames_to_column(var = "Gene_name") %>%
  dplyr::select(Gene_name, mean_score, source) %>%
  filter(mean_score > 0) %>%
  hgncConverter("Gene_name") %>%
  group_by(Gene_name) %>%
  mutate(mean_score = max(mean_score)) %>%
  distinct()


### Habermann et al. ----

habermann <- fread("files/banovich_markers.txt") %>%
  filter(
    cluster == "Ciliated" | cluster == "Differentiating Ciliated",
    p_val_adj < 0.001
  ) %>%
  mutate(
    mean_score = pct.1 - pct.2,
    Gene_name = gene
  ) %>%
  dplyr::select(Gene_name, mean_score) %>%
  group_by(Gene_name) %>%
  mutate(mean_score = exp(mean(log(mean_score)))) %>%
  filter(mean_score > 0) %>%
  hgncConverter("Gene_name") %>%
  group_by(Gene_name) %>%
  mutate(mean_score = max(mean_score)) %>%
  mutate(source = "habermann") %>%
  distinct()


### C. elegans single cell ----

celegans <- fread("files/celegans_sc_ciliated_oxygen_2022.txt") %>%
  filter(higher.expr == "Set 1" & q.val < 0.05 & (pct.1 - pct.2) > 0) %>%
  mutate(
    gene = as.character(gene),
    Gene_name = cele_orthology$Gene2Symbol[match(gene, cele_orthology$Gene1Symbol)],
    mean_score = normalization0(log(pct.1 - pct.2))
  ) %>%
  drop_na(Gene_name) %>%
  dplyr::select(15, 16) %>%
  group_by(Gene_name) %>%
  mutate(
    mean_score = exp(mean(log(mean_score))),
    source = "celegans"
  ) %>%
  group_by(Gene_name) %>%
  mutate(mean_score = max(mean_score)) %>%
  distinct()


### All Single cell 2 ----

all_sc<-rbind(carraro, habermann, reyfman, celegans) %>%
  dplyr::select(Gene_name) %>%
  distinct() %>% 
  mutate(carraro = ifelse(Gene_name %in% carraro$Gene_name, 1, 0),
         habermann = ifelse(Gene_name %in% habermann$Gene_name, 1, 0),
         reyfman = ifelse(Gene_name %in% reyfman$Gene_name, 1, 0),
         celegans = ifelse(Gene_name %in% celegans$Gene_name, 1, 0),
         ciliary = ifelse(Gene_name %in% ciliaryGenes$Gene_name, 1, 0),
         negative = ifelse(Gene_name %in% negativeCiliary$Gene_name, 1, 0)) %>%
  mutate(is.ciliary = ifelse(ciliary == 1, "Ciliary", ifelse(negative == 1, "Negative", "Unknown")))

all_sc2<-all_sc %>% 
  filter(ciliary == 1 | negative == 1)


# **********************
set.seed(444)
scindex<-sample(2, nrow(all_sc2), replace = TRUE, prob = c(0.7, 0.3))
scTrain<-all_sc2[scindex == 1,]
scTest<-all_sc2[scindex == 2,]
# **********************

library(ComplexUpset)
updata<-upset_data(all_sc2, colnames(all_sc)[2:5])
upset_test(all_sc2[,c(2:5,8)], colnames(all_sc2)[2:5])

sizes<-updata[['with_sizes']] %>%
  filter(in_exclusive_intersection == 1) %>%
  group_by(intersection) %>%
  mutate(count_ciliary = sum(is.ciliary == "Ciliary"),
         count_negative = sum(is.ciliary == "Negative")) %>%
  mutate(rate = count_ciliary / (count_ciliary + count_negative))

group_scores<-unique(sizes[,c(9,21)])


updata2<-upset_data(all_sc, colnames(all_sc)[2:5])
sizes2<-updata2[['with_sizes']] %>%
  filter(in_exclusive_intersection == 1) %>%
  left_join(group_scores, by = "intersection")
all_sc_scores<-sizes2[,c(1,19)]



sc2 <- rbind(carraro, reyfman, habermann, celegans) %>%
  group_by(Gene_name) %>%
  mutate(mean = mean(mean_score)) %>%
  distinct() %>%
  left_join(all_sc_scores, by = "Gene_name")
sc2$rate2<-impute.mean(sc2$rate)
sc2$score = sc2$mean*sc2$rate2

sc2$normalized_score<-normalization(sc2$score)
sc2<-unique(sc2[,c(1,8)])

xx<-left_join(sc2, ciliatest, by = "Gene_name")
xx$is.ciliary[is.na(xx$is.ciliary)]<-0

xxx<-left_join(single_cell_scores[,1], xx, by = "Gene_name")
xxx$normalized_score[is.na(xxx$normalized_score)]<-0
xxx$is.ciliary<-ifelse(xxx$Gene_name %in% negativeTest$Gene_name, 0, ifelse(xxx$Gene_name %in% ciliatest$Gene_name, 1, NA))
xxx<-xxx[!is.na(xxx$is.ciliary),]
classifier(xxx, "normalized_score", "is.ciliary")

library(pROC)
library(tidyverse)

rocs<-roc(is.ciliary ~ normalized_score, data = xxx)

# extract auc
rocs %>% 
  map(~tibble(AUC = .x$auc)) %>% 
  bind_rows(.id = "Methods") -> data.auc

# generate labels labels
data.auc %>% 
  mutate(label_long=paste0(Methods," , AUC = ",paste(round(AUC,2))),
         label_AUC=paste0("AUC = ",paste(round(AUC,2)))) -> data.labels

# plot on a single plot with AUC in labels
ggroc(rocs, size = 1.5) +
  scale_color_discrete(labels=data.labels$label_long) +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1),
               color="darkgrey", linetype="dashed")



### All Single cell ----

sc2 <- rbind(carraro, reyfman, habermann, celegans) %>%
  group_by(Gene_name) %>%
  mutate(ln = log(length(Gene_name) + 1)) %>%
  mutate(weighted_mean_score = ln * (exp(mean(log(mean_score + 1))) - 1)) %>%
  ungroup() %>%
  dplyr::select(Gene_name, weighted_mean_score) %>%
  distinct() %>%
  mutate(weighted_mean_score2 = normalization0(weighted_mean_score))


sc2 <- sc2 %>%
  mutate(score = weighted_mean_score2) %>%
  dplyr::select(Gene_name, score) %>%
  distinct()


single_cell_scores <- left_join(hgnc_names[, 2], sc2, by = c("Approved symbol" = "Gene_name"))
single_cell_scores$score[is.na(single_cell_scores$score)] <- 0
colnames(single_cell_scores)[1:2] <- c("Gene_name", "single_cell_score")


single_cell_scores$is.ciliary<-ifelse(single_cell_scores$Gene_name %in% negativeTest$Gene_name, 0, ifelse(single_cell_scores$Gene_name %in% ciliatest$Gene_name, 1, NA))
single_cell_scores<-single_cell_scores[!is.na(single_cell_scores$is.ciliary),]
classifier(single_cell_scores, "single_cell_score", "is.ciliary")

rocs<-roc(is.ciliary ~ single_cell_score, data = single_cell_scores)

write.table(single_cell_scores, "single_cell_scores.txt", sep = "\t", row.names = FALSE, quote = FALSE)
