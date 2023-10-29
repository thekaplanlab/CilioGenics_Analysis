

library(stringr)
library(biomaRt)
library(ontologyIndex)

### Biogrid ----

biogrid_mitab <- fread("files/BIOGRID-ALL-4.4.205.mitab.txt")
biogrid_mitab$`Confidence Values` <- as.numeric(gsub("score:", "", biogrid_mitab$`Confidence Values`))

# biogrid_mitab_pro<-biogrid_mitab[grepl("MI:0407|MI:0915|MI:0914|MI:0403", biogrid_mitab$`Interaction Types`),]
biogrid_mitab_pro <- biogrid_mitab

biogrid_mitab_pro_hum <- biogrid_mitab_pro %>%
  filter(`Taxid Interactor A` == "taxid:9606" & `Taxid Interactor B` == "taxid:9606") %>%
  mutate(
    Interactor_A = str_extract(`Alt IDs Interactor A`, "(?<=locuslink:).*?(?=\\||$)"),
    Interactor_B = str_extract(`Alt IDs Interactor B`, "(?<=locuslink:).*?(?=\\||$)")
  ) %>%
  dplyr::select(c(-1, -2, -3, -4, -5, -6)) %>%
  hgncConverter2("Interactor_B") %>%
  hgncConverter2("Interactor_A") %>%
  distinct() %>%
  mutate(Organism = "Homo sapiens") %>%
  relocate(Interactor_A, Interactor_B)


biogrid_mitab_pro_mus <- biogrid_mitab_pro %>%
  filter(`Taxid Interactor A` == "taxid:10090" & `Taxid Interactor B` == "taxid:10090") %>%
  mutate(
    Interactor_A = str_extract(`Alt IDs Interactor A`, "(?<=locuslink:).*?(?=\\|entrez|\\||$)"),
    Interactor_B = str_extract(`Alt IDs Interactor B`, "(?<=locuslink:).*?(?=\\|entrez|\\||$)")
  ) %>%
  dplyr::select(c(-1, -2, -3, -4, -5, -6)) %>%
  distinct()

biogrid_mitab_pro_mus$Interactor_A1 <- mouse_orthology2$`Approved symbol`[match(
  biogrid_mitab_pro_mus$Interactor_A,
  mouse_orthology2$Gene_name
)]

biogrid_mitab_pro_mus$Interactor_A2 <- mouse_orthology2$`Approved symbol`[match(
  biogrid_mitab_pro_mus$Interactor_A,
  mouse_orthology2$Gene_synonyms
)]

biogrid_mitab_pro_mus$Interactor_A <- ifelse(!is.na(biogrid_mitab_pro_mus$Interactor_A1),
                                             biogrid_mitab_pro_mus$Interactor_A1,
                                             biogrid_mitab_pro_mus$Interactor_A2
)

biogrid_mitab_pro_mus <- biogrid_mitab_pro_mus[, c(-12, -13)]

biogrid_mitab_pro_mus$Interactor_B1 <- mouse_orthology2$`Approved symbol`[match(
  biogrid_mitab_pro_mus$Interactor_B,
  mouse_orthology2$Gene_name
)]

biogrid_mitab_pro_mus$Interactor_B2 <- mouse_orthology2$`Approved symbol`[match(
  biogrid_mitab_pro_mus$Interactor_B,
  mouse_orthology2$Gene_synonyms
)]

biogrid_mitab_pro_mus$Interactor_B <- ifelse(!is.na(biogrid_mitab_pro_mus$Interactor_B1),
                                             biogrid_mitab_pro_mus$Interactor_B1,
                                             biogrid_mitab_pro_mus$Interactor_B2
)
biogrid_mitab_pro_mus <- biogrid_mitab_pro_mus[, c(-12, -13)]
biogrid_mitab_pro_mus$Organism <- "Mus musculus"



biogrid_mitab_pro_dro <- biogrid_mitab_pro %>%
  filter(`Taxid Interactor A` == "taxid:7227" & `Taxid Interactor B` == "taxid:7227") %>%
  mutate(
    Interactor_A = str_extract(`Alt IDs Interactor A`, "(?<=locuslink:).*?(?=\\|entrez|\\||$)"),
    Interactor_B = str_extract(`Alt IDs Interactor B`, "(?<=locuslink:).*?(?=\\|entrez|\\||$)")
  ) %>%
  dplyr::select(c(-1, -2, -3, -4, -5, -6)) %>%
  distinct()

biogrid_mitab_pro_dro$Interactor_A <- gsub('"', "", biogrid_mitab_pro_dro$Interactor_A)
biogrid_mitab_pro_dro$Interactor_B <- gsub('"', "", biogrid_mitab_pro_dro$Interactor_B)

biogrid_mitab_pro_dro$Interactor_A <- drosophila$Gene2Symbol[match(
  biogrid_mitab_pro_dro$Interactor_A,
  drosophila$Gene1Symbol
)]

biogrid_mitab_pro_dro$Interactor_B <- drosophila$Gene2Symbol[match(
  biogrid_mitab_pro_dro$Interactor_B,
  drosophila$Gene1Symbol
)]

biogrid_mitab_pro_dro$Organism <- "Drosophila melanogaster"

biogrid_mitab_pro_cele <- biogrid_mitab_pro %>%
  filter(`Taxid Interactor A` == "taxid:6239" & `Taxid Interactor B` == "taxid:6239") %>%
  mutate(
    Interactor_A = str_extract(`Alt IDs Interactor A`, "(?<=locuslink:).*?(?=\\|entrez|\\||$)"),
    Interactor_B = str_extract(`Alt IDs Interactor B`, "(?<=locuslink:).*?(?=\\|entrez|\\||$)")
  ) %>%
  dplyr::select(c(-1, -2, -3, -4, -5, -6)) %>%
  distinct()

biogrid_mitab_pro_cele$Interactor_A <- gsub("CELE_", "", biogrid_mitab_pro_cele$Interactor_A)
biogrid_mitab_pro_cele$Interactor_B <- gsub("CELE_", "", biogrid_mitab_pro_cele$Interactor_B)


biogrid_mitab_pro_cele$Interactor_A <- cele_orthology$Gene2Symbol[match(
  biogrid_mitab_pro_cele$Interactor_A,
  cele_orthology$Gene1Symbol
)]

biogrid_mitab_pro_cele$Interactor_B <- cele_orthology$Gene2Symbol[match(
  biogrid_mitab_pro_cele$Interactor_B,
  cele_orthology$Gene1Symbol
)]

biogrid_mitab_pro_cele$Organism <- "Caenorhabditis elegans"


mit_biogrid <- rbind(biogrid_mitab_pro_hum, biogrid_mitab_pro_mus, biogrid_mitab_pro_dro, biogrid_mitab_pro_cele)
mit_biogrid <- mit_biogrid[!is.na(mit_biogrid$Interactor_A) & !is.na(mit_biogrid$Interactor_B)]



### Intact ----

intact <- fread("files/intact.txt.1", select = c(5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 21, 22)) %>%
  unique() %>%
  mutate(ppi_score = as.numeric(str_extract(`Confidence value(s)`, "(?<=intact-miscore:).+"))) %>%
  filter(!is.na(ppi_score)) %>%
  mutate(`Publication Identifier(s)` = str_extract(`Publication Identifier(s)`, "pubmed:.*?(?=\\||$)"))


#### Intact human ----
intact_pro <- intact

intact_hum <- intact_pro %>% filter(
  `Taxid interactor A` == "taxid:9606(human)|taxid:9606(Homo sapiens)",
  `Taxid interactor B` == "taxid:9606(human)|taxid:9606(Homo sapiens)",
  grepl("gene name", `Alias(es) interactor A`),
  grepl("gene name", `Alias(es) interactor B`)
)

intact_hum1 <- intact_hum %>%
  mutate(`Alias(es) interactor A` = gsub(".*uniprotkb:(.+)\\(gene name).*", "\\1", intact_hum$`Alias(es) interactor A`)) %>%
  mutate(`Alias(es) interactor B` = gsub(".*uniprotkb:(.+)\\(gene name).*", "\\1", intact_hum$`Alias(es) interactor B`))

intact_hum1 <- hgncConverter2(intact_hum1, "Alias(es) interactor B")
intact_hum1 <- hgncConverter2(intact_hum1, "Alias(es) interactor A")

intact_hum1$Organism <- "Homo sapiens"

#### Intact celegans ----

intact_cele <- intact_pro %>% filter(
  `Taxid interactor A` == "taxid:6239(caeel)|taxid:6239(Caenorhabditis elegans)",
  `Taxid interactor B` == "taxid:6239(caeel)|taxid:6239(Caenorhabditis elegans)",
  grepl("orf name", `Alias(es) interactor A`) | grepl("gene name", `Alias(es) interactor A`),
  grepl("orf name", `Alias(es) interactor B`) | grepl("gene name", `Alias(es) interactor B`)
)

intact_cele[["Alias(es) interactor A"]] <- ifelse(grepl("\\(gene name\\)", intact_cele$`Alias(es) interactor A`),
                                                  gsub(".*uniprotkb:(.+)\\(gene name).*", "\\1", intact_cele$`Alias(es) interactor A`),
                                                  ifelse(grepl("CELE", intact_cele$`Alias(es) interactor A`),
                                                         str_extract(intact_cele$`Alias(es) interactor A`, "(?<=CELE_).+?(?=\\(orf name\\))"),
                                                         str_extract(intact_cele$`Alias(es) interactor A`, "(?<=uniprotkb:)(?!.*synonym).+?(?=\\(orf name\\))")
                                                  )
)


intact_cele[["Alias(es) interactor B"]] <- ifelse(grepl("\\(gene name\\)", intact_cele$`Alias(es) interactor B`),
                                                  gsub(".*uniprotkb:(.+)\\(gene name).*", "\\1", intact_cele$`Alias(es) interactor B`),
                                                  ifelse(grepl("CELE", intact_cele$`Alias(es) interactor B`),
                                                         str_extract(intact_cele$`Alias(es) interactor B`, "(?<=CELE_).+?(?=\\(orf name\\))"),
                                                         str_extract(intact_cele$`Alias(es) interactor B`, "(?<=uniprotkb:)(?!.*synonym).+?(?=\\(orf name\\))")
                                                  )
)

intact_cele1 <- intact_cele
intact_cele1$`Alias(es) interactor A` <- cele_orthology$Gene2Symbol[match(intact_cele1$`Alias(es) interactor A`, cele_orthology$Gene1Symbol)]
intact_cele1$`Alias(es) interactor B` <- cele_orthology$Gene2Symbol[match(intact_cele1$`Alias(es) interactor B`, cele_orthology$Gene1Symbol)]
intact_cele1 <- na.omit(intact_cele1)
intact_cele1$Organism <- "Caenorhabditis elegans"



#### Intact mouse ----

intact_mus <- intact %>% filter(
  `Taxid interactor A` == "taxid:10090(mouse)|taxid:10090(Mus musculus)",
  `Taxid interactor B` == "taxid:10090(mouse)|taxid:10090(Mus musculus)"
)

intact_mus1 <- intact_mus %>%
  mutate(`Alias(es) interactor A` = gsub(".*uniprotkb:(.+)\\(gene name).*", "\\1", intact_mus$`Alias(es) interactor A`)) %>%
  mutate(`Alias(es) interactor B` = gsub(".*uniprotkb:(.+)\\(gene name).*", "\\1", intact_mus$`Alias(es) interactor B`))

intact_mus1 <- intact_mus1 %>%
  mutate(`Alias(es) interactor A` = ifelse(nchar(intact_mus1$`Alias(es) interactor A`) < 15, `Alias(es) interactor A`,
                                           str_extract(`Alias(es) interactor A`, "(?<=psi-mi:).+?(?=_mouse_gene)")
  )) %>%
  mutate(`Alias(es) interactor B` = ifelse(nchar(intact_mus1$`Alias(es) interactor B`) < 15, `Alias(es) interactor B`,
                                           str_extract(`Alias(es) interactor B`, "(?<=psi-mi:).+?(?=_mouse_gene)")
  ))

intact_mus1 <- intact_mus1[!is.na(intact_mus1$`Alias(es) interactor A`) | !is.na(intact_mus1$`Alias(es) interactor A`), ]
intact_mus1$`Alias(es) interactor A` <- str_to_title(intact_mus1$`Alias(es) interactor A`)
intact_mus1$`Alias(es) interactor B` <- str_to_title(intact_mus1$`Alias(es) interactor B`)

intact_mus1$`Alias(es) interactor A` <- mouse_orthology2$`Approved symbol`[match(intact_mus1$`Alias(es) interactor A`, mouse_orthology2$Gene_name)]
intact_mus1$`Alias(es) interactor B` <- mouse_orthology2$`Approved symbol`[match(intact_mus1$`Alias(es) interactor B`, mouse_orthology2$Gene_name)]
intact_mus1 <- na.omit(intact_mus1)
intact_mus1$Organism <- "Mus musculus"


#### Intact drosophila ----

intact_dro <- intact %>%
  filter(`Taxid interactor A` == 'taxid:7227(drome)|taxid:7227("Drosophila melanogaster (Fruit fly)")') %>%
  filter(`Taxid interactor B` == 'taxid:7227(drome)|taxid:7227("Drosophila melanogaster (Fruit fly)")')

intact_dro[["Alias(es) interactor A"]] <- gsub(".*uniprotkb:(.+)\\(gene name).*", "\\1", intact_dro$`Alias(es) interactor A`)
intact_dro[["Alias(es) interactor A"]] <- gsub("Dmel\\\\", "", intact_dro$`Alias(es) interactor A`)
intact_dro[["Alias(es) interactor A"]] <- gsub("BEST:", "", intact_dro$`Alias(es) interactor A`)
intact_dro[["Alias(es) interactor A"]] <- gsub('"', "", intact_dro$`Alias(es) interactor A`)

intact_dro[["Alias(es) interactor B"]] <- gsub(".*uniprotkb:(.+)\\(gene name).*", "\\1", intact_dro$`Alias(es) interactor B`)
intact_dro[["Alias(es) interactor B"]] <- gsub("Dmel\\\\", "", intact_dro$`Alias(es) interactor B`)
intact_dro[["Alias(es) interactor B"]] <- gsub("BEST:", "", intact_dro$`Alias(es) interactor B`)
intact_dro[["Alias(es) interactor B"]] <- gsub('"', "", intact_dro$`Alias(es) interactor B`)


intact_dro1 <- intact_dro
intact_dro1$`Alias(es) interactor A` <- drosophila$Gene2Symbol[match(intact_dro1$`Alias(es) interactor A`, drosophila$Gene1Symbol)]
intact_dro1$`Alias(es) interactor B` <- drosophila$Gene2Symbol[match(intact_dro1$`Alias(es) interactor B`, drosophila$Gene1Symbol)]
intact_dro1 <- na.omit(intact_dro1)
intact_dro1$Organism <- "Drosophila melanogaster"



intact_all <- rbind(intact_hum1, intact_mus1, intact_dro1, intact_cele1)

onlyintact <- anti_join(intact_all, mit_biogrid, by = c(
  "Alias(es) interactor A" = "Interactor_A",
  "Alias(es) interactor B" = "Interactor_B",
  "Publication Identifier(s)" = "Publication Identifiers"
))

onlyintact <- onlyintact[, c(-10, -11, -12)]
mit_biogrid <- mit_biogrid[, -10]
colnames(onlyintact) <- colnames(mit_biogrid)

bioin <- rbind(mit_biogrid, onlyintact)


### Huri ----


huri <- fread("files/HI-union.psi", select = c(3, 4, 7, 8, 9, 10, 11, 12, 15), header = FALSE) %>%
  mutate(
    V3 = str_extract(V3, "ENSG.*?(?=\\.)"),
    V4 = str_extract(V4, "ENSG.*?(?=\\.)"),
    V15 = as.numeric(str_extract(V15, "(?<=author score:).*")),
    `Source Database` = "HuRI",
    Organism = "Homo sapiens"
  )

huri$V9[grepl("unassigned", huri$V9)] <- "pubmed:32296183"

# ensembl<-useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

hi_genes <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = unique(c(huri$V3, huri$V4)),
  mart = mart
) %>% na.omit()

hi_genes <- hgncConverter(hi_genes, "external_gene_name")
hi_genes$external_gene_name[hi_genes$external_gene_name == ""] <- hi_genes$ensembl_gene_id[hi_genes$external_gene_name == ""]

huri$V3 <- hi_genes$external_gene_name[match(huri$V3, hi_genes$ensembl_gene_id)]
huri$V4 <- hi_genes$external_gene_name[match(huri$V4, hi_genes$ensembl_gene_id)]

colnames(huri)[1:8] <- colnames(bioin)[1:8]
colnames(huri)[9] <- "Confidence Values"

onlyhuri <- anti_join(huri, bioin, by = c("Interactor_A", "Interactor_B", "Publication Identifiers"))

bioinhu <- rbind(bioin, onlyhuri)
bioinhu_pro <- bioinhu[!grepl("MI:0254", bioinhu$`Interaction Detection Method`), ]


### Scoring weights ----

owl <- get_ontology("files/mi.owl")

idmpar <- data.frame(
  method = c("MI:0013", "MI:0090", "MI:0254", "MI:0255", "MI:0401", "MI:0428"),
  weight = c(1, 0.66, 0.10, 0.10, 1, 0.33)
)

itpar <- data.frame(
  method = c("MI:0208", "MI:0403", "MI:0914", "MI:0915", "MI:0407"),
  weight = c(0.10, 0.33, 0.33, 0.66, 1),
  group = c(3, 1, 2, 2, 2)
)


### Genetic interactions ----

bioinhu_gen <- bioinhu[grepl("MI:0254", bioinhu$`Interaction Detection Method`), -10]

alliance_dro <- fread("files/dro_alliance.tsv.gz", select = c(5:13))
alliance_dro$Interactor_A <- str_extract(alliance_dro$`Alias(es) interactor A`, "(?<=flybase:).+?(?=\\(gene name\\))")
alliance_dro$Interactor_A <- gsub('"', "", alliance_dro$Interactor_A)
alliance_dro$Interactor_B <- str_extract(alliance_dro$`Alias(es) interactor B`, "(?<=flybase:).+?(?=\\(gene name\\))")
alliance_dro$Interactor_B <- gsub('"', "", alliance_dro$Interactor_B)
alliance_dro$Interactor_A <- drosophila$Gene2Symbol[match(alliance_dro$Interactor_A, drosophila$Gene1Symbol)]
alliance_dro$Interactor_B <- drosophila$Gene2Symbol[match(alliance_dro$Interactor_B, drosophila$Gene1Symbol)]
alliance_dro <- alliance_dro[!(is.na(alliance_dro$Interactor_A) | is.na(alliance_dro$Interactor_B)), ]

alliance_dro$`Publication Identifier(s)` <- str_extract(alliance_dro$`Publication Identifier(s)`, "pubmed:.+?(?=$)")
alliance_dro$`Publication Identifier(s)`[is.na(alliance_dro$`Publication Identifier(s)`)] <- "Unknown"

onlyalldro <- anti_join(alliance_dro, bioinhu_gen,
                        by = c(
                          "Publication Identifier(s)" = "Publication Identifiers",
                          "Interactor_A", "Interactor_B"
                        )
) %>%
  dplyr::select(c(10, 11, 3:9)) %>%
  mutate(Organism = "Drosophila melanogaster")

colnames(onlyalldro) <- colnames(bioinhu_gen)

alliance_cele <- fread("files/cele_alliance.tsv.gz", select = c(5:13))

alliance_cele$Interactor_A <- str_extract(alliance_cele$`Alias(es) interactor A`, "(?<=wormbase:).+?(?=\\(public_name\\))")
alliance_cele$Interactor_B <- str_extract(alliance_cele$`Alias(es) interactor B`, "(?<=wormbase:).+?(?=\\(public_name\\))")
alliance_cele$Interactor_A <- cele_orthology$Gene2Symbol[match(alliance_cele$Interactor_A, cele_orthology$Gene1Symbol)]
alliance_cele$Interactor_B <- cele_orthology$Gene2Symbol[match(alliance_cele$Interactor_B, cele_orthology$Gene1Symbol)]
alliance_cele <- alliance_cele[!(is.na(alliance_cele$Interactor_A) | is.na(alliance_cele$Interactor_B)), ]

onlyallcele <- anti_join(alliance_cele, bioinhu_gen,
                         by = c(
                           "Publication Identifier(s)" = "Publication Identifiers",
                           "Interactor_A", "Interactor_B"
                           )
                         ) %>%
  dplyr::select(c(10, 11, 3:9)) %>%
  mutate(Organism = "Caenorhabditis elegans")

colnames(onlyallcele) <- colnames(bioinhu_gen)

genetic_int_all <- rbind(bioinhu_gen, onlyalldro, onlyallcele)

bioinhu_pro <- rbind(bioinhu_pro[, -10], genetic_int_all)

bioinhu1 <- bioinhu_pro %>%
  mutate(
    method = str_extract(`Interaction Detection Method`, '(?<=psi-mi:").*?(?=")'),
    type = str_extract(`Interaction Types`, '(?<=psi-mi:").*?(?=")')
  )

### Score calculation ----

allmethods <- owl$ancestors[bioinhu1$method]
methodlist <- lapply(allmethods, function(x) idmpar$weight[which(idmpar$method %in% x)])
is.na(methodlist) <- lengths(methodlist) == 0
methodlist <- sapply(methodlist, "[[", 1)
methodlist <- unlist(methodlist, use.names = FALSE)

alltypes <- owl$ancestors[bioinhu1$type]
typelist <- lapply(alltypes, function(x) itpar$weight[which(itpar$method %in% x)])
is.na(typelist) <- lengths(typelist) == 0
typelist <- lapply(typelist, max)
typelist <- sapply(typelist, "[[", 1)
typelist <- unlist(typelist, use.names = FALSE)

grouplist <- lapply(alltypes, function(x) itpar$group[which(itpar$method %in% x)])
is.na(grouplist) <- lengths(grouplist) == 0
grouplist <- lapply(grouplist, max)
grouplist <- sapply(grouplist, "[[", 1)
grouplist <- unlist(grouplist, use.names = FALSE)


bioinhu1$idm <- methodlist
bioinhu1$it <- typelist
bioinhu1$it2 <- grouplist
bioinhu1$it[is.na(bioinhu1$it)] <- 0.05
bioinhu1$idm[is.na(bioinhu1$idm)] <- 0.05
bioinhu1$it2[is.na(bioinhu1$it2)] <- 3
bioinhu1 <- unique(bioinhu1)


yy1 <- bioinhu1 %>%
  group_by(Interactor_A, Interactor_B, `Interaction Detection Method`) %>%
  mutate(n = n_distinct(`Publication Identifiers`)) %>%
  ungroup() %>%
  group_by(Interactor_A, Interactor_B, `Interaction Types`) %>%
  mutate(nn = n_distinct(`Publication Identifiers`)) %>%
  ungroup() %>%
  mutate(
    a = idm * n,
    b = it * nn
  )

yy1 <- yy1[!is.na(yy1$Interactor_A), ]
yy1 <- yy1[!is.na(yy1$Interactor_B), ]

yy2 <- yy1 %>%
  group_by(Interactor_A, Interactor_B) %>%
  distinct(`Interaction Detection Method`, .keep_all = TRUE) %>%
  mutate(a1 = sum(a))

yy3 <- yy1 %>%
  group_by(Interactor_A, Interactor_B) %>%
  distinct(`Interaction Types`, .keep_all = TRUE) %>%
  mutate(b1 = sum(b))

yy4 <- yy1 %>%
  left_join(yy2, by = names(yy1)) %>%
  left_join(yy3, by = names(yy1)) %>%
  group_by(Interactor_A, Interactor_B) %>%
  mutate(a1 = max(a1, na.rm = TRUE)) %>%
  mutate(b1 = max(b1, na.rm = TRUE)) %>%
  group_by(Interactor_A, Interactor_B, it2) %>%
  mutate(gscv2 = sum(max(it))) %>%
  group_by(Interactor_A, Interactor_B) %>%
  mutate(gscv = sum(max(idm))) %>%
  ungroup() %>%
  mutate(method_score = log(a1 + 1, a1 + gscv + 1)) %>%
  mutate(type_score = log(b1 + 1, b1 + gscv2 + 1))

yy5 <- yy4[, c(-11:-23)]
# yy5$ppi_score<-(sqrt(yy5$method_score * yy5$type_score))
yy5$ppi_score <- (yy5$method_score + yy5$type_score) / 2
yy5 <- as.data.frame(yy5)

# **************************************

# yy5$c <- ifelse(yy5$Interactor_B %in% ciliaryGenes$Gene_name, 1, 0)
# yy5$c2 <- ifelse(yy5$Interactor_A %in% ciliaryGenes$Gene_name, 1, 0)
# 
# yy6 <- yy5 %>%
#   group_by(Interactor_A) %>%
#   mutate(
#     rx = sum(ppi_score),
#     rcx = sum(ppi_score[c == 1])
#   ) %>%
#   mutate(sout = rcx / rx) %>%
#   ungroup() %>%
#   group_by(Interactor_B) %>%
#   mutate(
#     px = sum(1 / rx),
#     pcx = sum(1 / rx[c2 == 1])
#   ) %>%
#   mutate(sin = pcx / px)
# 
# yy61 <- yy6[, c(1, 18)]
# yy62 <- yy6[, c(2, 21)]
# 
# colnames(yy61) <- c("Gene_name", "scores")
# colnames(yy62) <- c("Gene_name", "scores")
# 
# all_interactions <- rbind(yy61, yy62) %>%
#   group_by(Gene_name) %>%
#   mutate(interaction_scores = sum(scores)) %>%
#   ungroup() %>%
#   dplyr::select(Gene_name, interaction_scores) %>%
#   distinct() %>%
#   mutate(sc = log(interaction_scores + 1) - 1)
# 
# all_interactions$sc <- normalization0(all_interactions$sc)
# colnames(all_interactions)[3] <- "interaction_scores"
# all_interactions <- all_interactions[, c(1, 3)]
# prot_int1 <- left_join(hgnc_names[, 2], all_interactions, by = c("Approved symbol" = "Gene_name"))
# colnames(prot_int1) <- c("Gene_name", "interaction_scores")
# prot_int$interaction_scores[is.na(prot_int$interaction_scores)] <- 0
# 
# write.table(prot_int, "protein_interaction_scores_2022v3_o.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# **************************************

yy6 <- yy5 %>%
  left_join(ciliaryGenes, by = c("Interactor_B" = "Gene_name")) %>%
  mutate(across(c("is.ciliary"), ~ replace_na(.x, 0))) %>%
  rename("Ciliary_A" = "is.ciliary")

yy6 <- yy6 %>%
  left_join(ciliaryGenes, by = c("Interactor_A" = "Gene_name")) %>%
  mutate(across(c("is.ciliary"), ~ replace_na(.x, 0))) %>%
  rename("Ciliary_B" = "is.ciliary") %>%
  mutate(
    Ciliary_A = Ciliary_A * ppi_score,
    Ciliary_B = Ciliary_B * ppi_score
  )


int_a <- unique(yy6[, c(1, 2, 14)])
int_b <- unique(yy6[, c(2, 1, 15)])
colnames(int_b) <- colnames(int_a)

int_all <- rbind(int_a, int_b)

alist <- int_all[, c(1, 3)]
alist1x <- alist %>%
  group_by(Interactor_A) %>%
  mutate(gscv = sum(max(Ciliary_A)), a1 = sum(Ciliary_A), ln = length(Interactor_A)) %>%
  mutate(ln1 = log((ln / 5) + 1)) %>%
  mutate(scoreFinal = ln1 * (exp(mean(log(Ciliary_A + 1))) - 1))


alist2x <- unique(alist1x[, c(1, 7)])
colnames(alist2x) <- c("Gene_name", "interaction_scores")
alist2x$interaction_scores <- log(alist2x$interaction_scores + 1) - 1
alist2x$interaction_scores <- normalization0(alist2x$interaction_scores)

alist2x$interaction_scores[alist2x$interaction_scores != 0] <-
  normalization0(log(alist2x$interaction_scores[alist2x$interaction_scores != 0]))

prot_int <- left_join(hgnc_names[, 2], alist2x, by = c("Approved symbol" = "Gene_name"))
colnames(prot_int) <- c("Gene_name", "interaction_scores")
prot_int$interaction_scores[is.na(prot_int$interaction_scores)] <- 0

write.table(prot_int, "protein_interaction_scores", row.names = FALSE, quote = FALSE, sep = "\t")
