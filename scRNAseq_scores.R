
source("data.R")

## Single cell scores ----

### Reyfman ----

reyfman <- readRDS("./files/markers_reyfman.RDS") %>%
  filter(
    cluster == "Ciliated Cells", pct.1 - pct.2 > 0,
    p_val_adj < 0.001
  ) %>%
  hgncConverter("gene") %>%
  mutate(
    #mean_score = normalization0(pct.1 - pct.2),
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

carraro <- readRDS("./files/markers.RDS")
carraro <- carraro[[2]] %>%
  filter(p_val_adj < 0.001, pct.1 - pct.2 > 0) %>%
  mutate(
    mean_score = normalization0(pct.1 - pct.2),
    source = "carraro"
  ) %>%
  tibble::rownames_to_column(var = "Gene_name") %>%
  dplyr::select(Gene_name, mean_score, source) %>%
  hgncConverter("Gene_name") %>%
  group_by(Gene_name) %>%
  mutate(mean_score = max(mean_score)) %>%
  distinct()


### Habermann et al. ----

habermann <- fread("./files/banovich_markers.txt") %>%
  filter(
    cluster == "Ciliated" | cluster == "Differentiating Ciliated",
    p_val_adj < 0.001, pct.1 - pct.2 > 0
  ) %>%
  mutate(
    mean_score = normalization0(pct.1 - pct.2),
    Gene_name = gene
  ) %>%
  dplyr::select(Gene_name, mean_score) %>%
  group_by(Gene_name) %>%
  mutate(mean_score = exp(mean(log(mean_score)))) %>%
  hgncConverter("Gene_name") %>%
  group_by(Gene_name) %>%
  mutate(mean_score = max(mean_score)) %>%
  mutate(source = "habermann") %>%
  distinct()


# Murthy et al.

murthy <- readRDS("./files/murthy_diff.RDS")
#murthy$Gene_name<-rownames(murthy)
murthy <- murthy %>%
  filter(p_val_adj < 0.001, pct.1 - pct.2 > 0) %>%
  mutate(
    mean_score = normalization0(pct.1 - pct.2),
    #mean_score = normalization0(log(pct.1 - pct.2)),
    #mean_score = pct.1 - pct.2,
    source = "murthy"
  ) %>%
  tibble::rownames_to_column(var = "Gene_name") %>%
  dplyr::select(Gene_name, mean_score, source) %>%
  #filter(mean_score > 0.2) %>%
  hgncConverter("Gene_name") %>%
  distinct()


### C. elegans single cell ----

celegans <- fread("./files/celegans_sc_ciliated_oxygen_2022.txt") %>%
  filter(higher.expr == "Set 1" & q.val < 0.05 & (pct.1 - pct.2) > 0) %>%
  mutate(
    gene = as.character(gene),
    Gene_name = cele_orthology$Gene2Symbol[match(gene, cele_orthology$Gene1Symbol)],
    #mean_score = normalization0(pct.1 - pct.2)
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


### Trachea single cell ----

trachea <- fread("./files/trachea_all_markers.txt") %>%
  filter(
    cluster == "Ciliated Cells", pct.1 - pct.2 > 0,
    p_val_adj < 0.001
  ) %>%
  hgncConverter("gene_name") %>%
  mutate(
    mean_score = normalization0(pct.1 - pct.2),
    source = "trachea"
  ) %>%
  rename("Gene_name" = "gene_name") %>%
  dplyr::select(Gene_name, mean_score, source) %>%
  group_by(Gene_name) %>%
  mutate(mean_score = max(mean_score)) %>%
  distinct()
trachea$mean_score[trachea$mean_score == 0] <- 0.00001
trachea<-trachea[-c(1,2),]


### CenGen single cell ----

cengen <- fread("./files/difex_cengen.csv") %>%
  filter(human_symbol != "") %>%
  filter(
    cluster == "Sensory", pct.1 - pct.2 > 0,
    p_val_adj < 0.001
  ) %>%
  hgncConverter("human_symbol") %>%
  mutate(
    mean_score = normalization0(log(pct.1 - pct.2)),
    source = "cengen"
  ) %>%
  rename("Gene_name" = "human_symbol") %>%
  dplyr::select(Gene_name, mean_score, source) %>%
  group_by(Gene_name) %>%
  mutate(mean_score = max(mean_score)) %>%
  distinct()
cengen$mean_score[cengen$mean_score == 0] <- 0.00001


### All Single cell ----

sc2 <- rbind(carraro, reyfman, habermann, celegans, murthy, trachea, cengen) %>%
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

write.table(single_cell_scores, "single_cell_scores.txt", sep = "\t", row.names = FALSE, quote = FALSE)
