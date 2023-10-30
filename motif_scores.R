
source("data.R")

## Motif scores ----

pan_motif <- read_xlsx("./files/Table S16.xlsx", sheet = "Ciliated_cells") %>%
  dplyr::select(total, pairs)
pan_motif$motif <- str_split(pan_motif$pairs, "_", simplify = TRUE)[, 1]
pan_motif$Gene_name <- str_split(pan_motif$pairs, "_", simplify = TRUE)[, 2]

pan_motif <- left_join(pan_motif, mouse_orthology2[, c(2, 7)], by = c("motif" = "Gene_name"))
pan_motif$Gene_name2 <- mouse_orthology2$`Approved symbol`[match(pan_motif$Gene_name, mouse_orthology2$Gene_name)]
pan_motif$Gene_name3 <- mouse_orthology2$`Approved symbol`[match(pan_motif$Gene_name, mouse_orthology2$Gene_synonyms)]

pan_motif$Gene_name <- ifelse(!is.na(pan_motif$Gene_name2), pan_motif$Gene_name2, pan_motif$Gene_name3)
pan_motif <- pan_motif[!is.na(pan_motif$Gene_name), c(1, 4, 5)]
colnames(pan_motif)[3] <- "motif"
pan_motif <- unique(pan_motif)
pan_motif <- hgncConverter(pan_motif, "Gene_name")


## foxj1 ***********************

foxj1 <- read_xlsx("./files/SupplemantaryTable1.xlsx", sheet = 2) %>%
  dplyr::select(`Target Gene Name`) %>%
  count(`Target Gene Name`) %>%
  rename("Gene_name" = `Target Gene Name`, "total" = n) %>%
  mutate(motif = "FOXJ1") %>%
  hgncConverter("Gene_name")

# foxj1$total<-(normalization0(foxj1$total)*6)+1

pan_motif <- rbind(pan_motif, foxj1)

motif_list<-c("RFX2", "RFX3", "MYB", "GLIS3", "JAZF1", "FOXJ1", "SOX5", "TOX")

cil_motif<-pan_motif[pan_motif$Gene_name %in% ciliaryGenes$Gene_name,]

pan_motifx<-pan_motif[pan_motif$motif %in% motif_list,]

motif_scores<-pan_motifx %>% 
  mutate(total2 = log(total + 1)) %>%
  dplyr::select(Gene_name, total2) %>%
  group_by(Gene_name) %>%
  mutate(motif_score = sum(total2)) %>%
  dplyr::select(Gene_name, motif_score) %>%
  distinct() %>%
  mutate(motif_score = log(motif_score + 1))

motif_scores2 <- motif_scores
motif_scores2 <- left_join(hgnc_names[, 2], motif_scores2, by = c("Approved symbol" = "Gene_name"))
colnames(motif_scores2)[1]<-"Gene_name"
motif_scores2[is.na(motif_scores2),]<-0

motif_scores2$motif_score<-normalization(motif_scores2$motif_score)

write.table(motif_scores2, "motif_scores.txt", sep = "\t", row.names = FALSE, quote = FALSE)

