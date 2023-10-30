


source("interaction_scores.R")
source("scRNAseq_scores.R")
source("comparative_genomics_scores.R")
source("motif_scores.R")
source("protein_atlas_scores.R")

## CilioGenics score calculation ----

interaction_scores <- fread("protein_interaction_scores.txt")
cluster_scores <- fread("cluster_scores.txt")
motif_scores <- fread("motif_scores.txt")
proatlas_scores <- fread("proatlas_score.txt")
sc_scores <- fread("single_cell_scores.txt")


all_scores <- left_join(hgnc_names[, 2], interaction_scores, by = c("Approved symbol" = "Gene_name"))
colnames(all_scores)[1] <- "Gene_name"
colnames(cluster_scores)[1] <- "Gene_name"
all_scores <- left_join(all_scores, cluster_scores, by = "Gene_name")
all_scores <- left_join(all_scores, motif_scores, by = "Gene_name")
all_scores <- left_join(all_scores, proatlas_scores, by = "Gene_name")
all_scores <- left_join(all_scores, sc_scores, by = "Gene_name")

