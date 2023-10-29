

## Motif scores ----

pan_motif <- read_xlsx("files/Table S16.xlsx", sheet = "Ciliated_cells") %>%
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

foxj1 <- read_xlsx("files/SupplemantaryTable1.xlsx", sheet = 2) %>%
  dplyr::select(`Target Gene Name`) %>%
  count(`Target Gene Name`) %>%
  rename("Gene_name" = `Target Gene Name`, "total" = n) %>%
  mutate(motif = "FOXJ1") %>%
  hgncConverter("Gene_name")

# foxj1$total<-(normalization0(foxj1$total)*6)+1

pan_motif <- rbind(pan_motif, foxj1)

motif_list<-c("RFX2", "RFX3", "MYB", "GLIS3", "JAZF1", "FOXJ1", "SOX5", "TOX")
non_motif_list<-c("FOXO3", "CUX1", "CREB3", "PRDM6", "ETV5", "MITF", "ZEB2", "BNC2", "NR2C1")

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

motif_scores3<-motif_scores2[motif_scores2$Gene_name %in% c(ciliaryGenes$Gene_name, negativeCiliary$Gene_name),]
motif_scores3$is.ciliary<-ifelse(motif_scores3$Gene_name %in% ciliaryGenes$Gene_name, 1, 0)

classifier(motif_scores3, "motif_score", "is.ciliary")

motif_scores3<-motif_scores2[motif_scores2$Gene_name %in% c(ciliatest$Gene_name, negativeCiliary$Gene_name),]
motif_scores3$is.ciliary<-ifelse(motif_scores3$Gene_name %in% ciliatest$Gene_name, 1, 0)

motif_scores3<-motif_scores2[motif_scores2$Gene_name %in% c(b, negativeCiliary$Gene_name),]
motif_scores3$is.ciliary<-ifelse(motif_scores3$Gene_name %in% b, 1, 0)




# figure #













tfs <- data.frame(
  motif = c("RFX2", "RFX3", "MYB", "GLIS3", "JAZF1", "FOXJ1", "ZBTB16", "SOX5", "TOX"),
  score = c(1, 1, 1, 1, 1, 1, 1, 1, 1)
)



# ******************************

pan_motif$predicted <- ifelse(pan_motif$motif %in% c("RFX2", "RFX3", "MYB", "GLIS3", "JAZF1", "FOXJ1", "ZBTB16", "SOX5", "TOX"), 1, 0)
pan_motif$is.ciliary <- ifelse(pan_motif$Gene_name %in% ciliaryGenes$Gene_name, 1, 0)

pan_test <- pan_motif %>%
  mutate(
    predicted = as.factor(predicted),
    is.ciliary = as.factor(is.ciliary)
  )

### ***************************
pan_motif$is.ciliary <- ifelse(pan_motif$Gene_name %in% ciliaryGenes$Gene_name, 1,
                               ifelse(pan_motif$Gene_name %in% negativeTrain$Gene_name, 0, NA)
)
pan_train <- pan_motif[!is.na(pan_motif$is.ciliary), ]
pan_test <- pan_train %>%
  mutate(
    predicted = as.factor(predicted),
    is.ciliary = as.factor(is.ciliary)
  )

# tfs <- data.frame(
#   motif = c("RFX2", "RFX3", "MYB", "GLIS3", "JAZF1", "FOXJ1"),
#   score = c(1, 0.965, 1, 0.970, 0.973, 0.973)
# )

### ***************************

cnfs <- confusionMatrix(pan_test$predicted[pan_test$motif %in% c("RFX2", "RFX3", "MYB", "GLIS3", "JAZF1", "FOXJ1", "ZBTB16", "SOX5", "TOX")],
                        pan_test$is.ciliary[pan_test$motif %in% c("RFX2", "RFX3", "MYB", "GLIS3", "JAZF1", "FOXJ1", "ZBTB16", "SOX5", "TOX")],
                        positive = "1", mode = "prec_recall"
)
# cnfs <- confusionMatrix(pan_test$predicted[pan_test$motif == "JAZF1"],
#                         pan_test$is.ciliary[pan_test$motif == "JAZF1"],
#                         positive = "1", mode = "prec_recall"
# )
paste0("FOXJ1 is ", cnfs$byClass[["Precision"]])
# classifier(pan_motif[pan_motif$motif == "RFX2",], "predicted", "is.ciliary")

pan_rate <- length(unique(pan_motif$Gene_name[pan_motif$is.ciliary == 1 & pan_motif$predicted == 1])) / length(ciliaryGenes$Gene_name)
pan_rate2 <- length(pan_motif$Gene_name[pan_motif$is.ciliary == 1 & pan_motif$predicted == 1]) / length(pan_motif$Gene_name[pan_motif$predicted == 1])
pan_motif$score <- pan_motif$total * pan_motif$predicted

# tfs<-data.frame(motif = c("RFX2","RFX3","MYB","GLIS3","JAZF1"),
#                 score = c(0.15,0.26,0.08,0.19,0.14))

tfs <- data.frame(
  motif = c("RFX2", "RFX3", "MYB", "GLIS3", "JAZF1"),
  score = c(0.576, 0.336, 0.568, 0.369, 0.497)
)

pan_motif$scv <- tfs$score[match(pan_motif$motif, tfs$motif)]
pan_motif$scv[is.na(pan_motif$scv)] <- 0

pan_motif1 <- pan_motif %>%
  group_by(Gene_name) %>%
  mutate(a = sum(scv * total)) %>%
  mutate(b = a + sum(max(scv))) %>%
  mutate(score = log(a + 1, b + 1)) %>%
  dplyr::select(Gene_name, score) %>%
  distinct()
pan_motif1$score[is.nan(pan_motif1$score)] <- 0 # V3
pan_motif1[nrow(pan_motif1) + 1, 2] <- min(pan_motif1$score[!is.nan(pan_motif1$score)]) - 0.1

# pan_motif1$score[is.nan(pan_motif1$score)]<-0

# pan_motif2<-pan_motif[pan_motif$motif %in% c("RFX2","RFX3","MYB","GLIS3","JAZF1"),]
# pan_motif2<-pan_motif %>%
#   mutate(a = scv * log(total)) %>%
#   group_by(Gene_name) %>%
#   mutate(b = (log(length(Gene_name)+1)-1)*(exp(mean(log(a))))) %>%
#   mutate(score = b) %>%
#   dplyr::select(Gene_name, score) %>%
#   distinct()

pan_motif1$score <- normalization0(pan_motif1$score)
pan_motif1$score[is.nan(pan_motif1$score)] <- 0
pan_motif1 <- pan_motif1[!is.na(pan_motif1$Gene_name), ]

# pan_motif2$score<-normalization0(pan_motif2$score)
pan_motif2 <- pan_motif1
pan_motif2 <- left_join(hgnc_names[, 2], pan_motif2, by = c("Approved symbol" = "Gene_name"))
colnames(pan_motif2) <- c("Gene_name", "motif_score")
pan_motif2$motif_score[is.na(pan_motif2$motif_score)] <- 0
motif_scores <- pan_motif2

write.table(motif_scores, "motif_scores.txt", sep = "\t", row.names = FALSE, quote = FALSE)
