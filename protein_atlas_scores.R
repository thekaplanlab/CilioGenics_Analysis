



## Protein atlas scores ----

#source("proteinAtlas_webcrawling.R")

proatlas2 <- fread("./files/protein_atlas_2022.txt")
colnames(proatlas2) <- c(
  "Gene_name", "Cilia_comments", "Cilium_comments", "Centrosome_comments", "Flagella_comments", "Flagellum_comments",
  "Cilia", "Cilium", "Centrosome", "Flagella", "Flagellum"
)

proatlas3 <- subset(proatlas2, Cilia == "YES" | Cilium == "YES" | Centrosome == "YES" | Flagella == "YES" | Flagellum == "YES")

progenes <- proatlas3[, 1]

progenes$Centrosome <- ifelse(apply(proatlas3[, 2:6], 1, function(x) any(grepl("centrosome", x, ignore.case = TRUE))), "YES", "NO")
progenes$Cilia_expression <- ifelse(apply(proatlas3[, 2:6], 1, function(x) any(grepl("cilia|cilium", x, ignore.case = TRUE))), "YES", "NO")
progenes$Flagella <- ifelse(apply(proatlas3[, 2:6], 1, function(x) any(grepl("flagella|flagellum", x, ignore.case = TRUE))), "YES", "NO")
progenes$Cilia_localization <- ifelse(apply(proatlas3[, 2:6], 1, function(x) any(grepl("cilia|cilium", x, ignore.case = TRUE) & grepl("positivity|stain", x, ignore.case = TRUE))), "YES", "NO")

progenes$Centrosome <- ifelse(progenes$Centrosome == "YES", 0.25, 0)
progenes$Cilia_expression <- ifelse(progenes$Cilia_expression == "YES", 0.5, 0)
progenes$Flagella <- ifelse(progenes$Flagella == "YES", 1, 0)
progenes$Cilia_localization <- ifelse(progenes$Cilia_localization == "YES", 1, 0)

progenes$score <- apply(progenes[, 2:5], 1, max)
progenes$score <- apply(progenes[, 2:5], 1, sum)
progenes <- hgncConverter(progenes, "Gene_name") # This
progenes<-progenes %>%
  group_by(Gene_name) %>%
  summarise_all(sum)


protein_atlas_score <- progenes[, c(1, 6)]

protein_atlas_score <- left_join(hgnc_names[, 2], protein_atlas_score, by = c("Approved symbol" = "Gene_name"))
protein_atlas_score$score[is.na(protein_atlas_score$score)] <- 0
protein_atlas_score$score<-normalization(protein_atlas_score$score)
colnames(protein_atlas_score)[1:2] <- c("Gene_name", "pa_score")

write.table(protein_atlas_score, "proatlas_score.txt", sep = "\t", row.names = FALSE, quote = FALSE)

