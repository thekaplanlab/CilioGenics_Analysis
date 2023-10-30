
# All ----

## Generate Objects ----

## Load libraries ----

library(data.table)
library(geneName)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(biomaRt)
library(ontologyIndex)
library(cluster)


yannis<-read_xlsx(paste0("files/SupplementaryTable2_edited.xlsx")) %>%
  rename("Gene_name" = "Gene Name") %>%
  mutate(yannis = 1)
clime<-fread(paste0("files/human_cilia_clime.txt"), select = "Gene Symbol") %>%
  rename("Gene_name" = "Gene Symbol") %>%
  mutate(clime = 1)


### Load functions ----

source("functions.R")

### Read files required for common ----

hgnc_names<-fread("files/hgnc_names.txt", select = c(1,3,11)) %>%
  filter(!grepl("-", `Approved symbol`)) %>%
  dplyr::select(1,2) %>%
  unique()

# hgnc_names<-fread("hgnc_complete_set.txt") %>%
#   filter(locus_group == "protein-coding gene",
#          status == "Approved") %>%
#   dplyr::select(hgnc_id, symbol) %>%
  

### Ciliary and non ciliary gene lists

ciliaryGenes<-fread("files/goldstandart_v2.txt")
ciliaryGenes$is.ciliary<-1
ciliaryGenes<-hgncConverter(ciliaryGenes, "Gene_name")
ciliaryGenes$Gene_name[ciliaryGenes$Gene_name == "CILK"]<-"CILK1"
ciliaryGenes$Gene_name[ciliaryGenes$Gene_name == "GPBAR"]<-"GPBAR1"
ciliaryGenes<-na.omit(ciliaryGenes)

negativeCiliary<-read_xls("files/Nevers_2017_NegativeGenesInsightsCiiaryGenes_SuppTable3.xls")
negativeCiliary<-negativeCiliary[,1]
colnames(negativeCiliary)<-"Gene_name"
negativeCiliary<-hgncConverter(negativeCiliary, "Gene_name")
negativeCiliary$is.ciliary<-0

negativeCiliary<-anti_join(negativeCiliary, ciliaryGenes, by = "Gene_name")

ciliatest<-read_xlsx("files/cilia_genes_references.xlsx")
colnames(ciliatest)[1]<-"Gene_name"
ciliatest<-hgncConverter(ciliatest, "Gene_name")
ciliatest<-anti_join(ciliatest, ciliaryGenes, by = "Gene_name")

negativeCiliary<-anti_join(negativeCiliary, ciliatest, by = "Gene_name")
# set.seed(444)
# negativeIndex<-sample(2, nrow(negativeCiliary), replace = TRUE, prob = c(0.5, 0.5))
# negativeTrain<-negativeCiliary[negativeIndex == 1,]
# negativeTest<-negativeCiliary[negativeIndex == 2,]

# **********************************
set.seed(4444)
negativeTest<-data.frame(Gene_name = sample(negativeCiliary$Gene_name, 487), is.ciliary = 0)
negativeTrain<-anti_join(negativeCiliary, negativeTest, by = "Gene_name")

# **********************************


ciliatest<-ciliatest %>% dplyr::select(Gene_name) %>% mutate(is.ciliary = 1)

alltest<-rbind(ciliatest, negativeTest)


### Orthology data ----

orthology<-fread("files/ORTHOLOGY-ALLIANCE_COMBINED_0.tsv.gz") %>% filter(IsBestScore == "Yes")
orthology<-orthology[orthology$Gene2SpeciesName == "Homo sapiens", c(2,4,6,8)]
orthology<-hgncConverter(orthology, "Gene2Symbol")

#### Mouse ----

mouse_orthology<-fread("files/HGNC_AllianceHomology.rpt")
colnames(mouse_orthology)[1:16]<-colnames(mouse_orthology)[2:17]
mouse_orthology<-mouse_orthology[,-17]
colnames(mouse_orthology)[2]<-"Gene_name"

mouse_synonyms<-fread("files/MRK_List2.rpt")
colnames(mouse_synonyms)[c(7,12)]<-c("Gene_name","Gene_synonyms")
mouse_synonyms<-separate_rows(mouse_synonyms, Gene_synonyms, sep = "\\|") %>% 
  dplyr::select(7,12) %>% unique()
mouse_synonyms<-mouse_synonyms[mouse_synonyms$Gene_synonyms != "",]
mouse_orthology<-left_join(mouse_orthology, mouse_synonyms, by = "Gene_name")

hgnc_mouse<-fread("files/hgnc_names.txt") %>% 
  filter(Status == "Approved") %>% 
  dplyr::select(1,3) %>% 
  hgncConverter("Approved symbol")

mouse_orthology2<-separate_rows(mouse_orthology, `HGNC ID`, sep = "\\|")
mouse_orthology2<-mouse_orthology2[,c(1,2,5,10,16,17)]
mouse_orthology2<-left_join(mouse_orthology2, hgnc_mouse, by = "HGNC ID")


#### C. elegans ----

cele_orthology<-orthology %>% filter(Gene1SpeciesName == "Caenorhabditis elegans")
cele_ids<-fread("https://downloads.wormbase.org/releases/WS285/species/c_elegans/PRJNA13758/annotation/c_elegans.PRJNA13758.WS285.geneIDs.txt.gz") %>% 
  mutate(V3 = ifelse(V3 == "", V4, V3)) %>%
  rename("Gene1Symbol" = "V3", "Gene_ID" = "V2") %>% 
  dplyr::select(Gene1Symbol, Gene_ID) %>% 
  left_join(cele_orthology[,c(1,2)], by = "Gene1Symbol")

#### Drosophila ----

drosophila<-fread("files/dmel_human_orthologs_disease_fb_2021_06.tsv.gz")
ort_dro<-orthology[orthology$Gene1SpeciesName == "Drosophila melanogaster",]
ort_dro$Gene1Symbol<-gsub("α","alpha",ort_dro$Gene1Symbol)
onlyortdro<-anti_join(ort_dro, drosophila, by = c("Gene1Symbol" = "Dmel_gene_symbol"))

drosophila<-drosophila[,c(2,5)]
colnames(drosophila)<-c("Gene1Symbol","Gene2Symbol")
drosophila<-rbind(drosophila, onlyortdro[,c(1,2)]) %>% unique() %>% hgncConverter("Gene2Symbol")


## Interactions ----

### Biogrid ----

biogrid_mitab <- fread("files/BIOGRID-ALL-4.4.205.mitab.txt")
biogrid_mitab$`Confidence Values` <- as.numeric(gsub("score:", "", biogrid_mitab$`Confidence Values`))

# biogrid_mitab_pro<-biogrid_mitab[grepl("MI:0407|MI:0915|MI:0914|MI:0403", biogrid_mitab$`Interaction Types`),]
biogrid_mitab_pro <- as.data.frame(biogrid_mitab)

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
mit_biogrid <- mit_biogrid[!is.na(mit_biogrid$Interactor_A) & !is.na(mit_biogrid$Interactor_B),]



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

write.table(prot_int, "protein_interaction_scores.txt", row.names = FALSE, quote = FALSE, sep = "\t")


## Single cell scores ----

### Reyfman ----

reyfman <- readRDS("files/markers_reyfman.RDS") %>%
  filter(
    cluster == "Ciliated Cells", pct.1 - pct.2 > 0,
    p_val_adj < 0.001
  ) %>%
  hgncConverter("gene") %>%
  mutate(
    mean_score = normalization0(pct.1 - pct.2),
    #mean_score = normalization0(log(pct.1 - pct.2)),
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
  filter(p_val_adj < 0.001, pct.1 - pct.2 > 0) %>%
  mutate(
    mean_score = normalization0(pct.1 - pct.2),
    #mean_score = normalization0(log(pct.1 - pct.2)),
    #mean_score = pct.1 - pct.2,
    source = "carraro"
  ) %>%
  tibble::rownames_to_column(var = "Gene_name") %>%
  dplyr::select(Gene_name, mean_score, source) %>%
  hgncConverter("Gene_name") %>%
  group_by(Gene_name) %>%
  mutate(mean_score = max(mean_score)) %>%
  distinct()


### Habermann et al. ----

habermann <- fread("files/banovich_markers.txt") %>%
  filter(
    cluster == "Ciliated" | cluster == "Differentiating Ciliated",
    p_val_adj < 0.001, pct.1 - pct.2 > 0
  ) %>%
  mutate(
    mean_score = normalization0(pct.1 - pct.2),
    #mean_score = normalization0(log(pct.1 - pct.2)),
    #mean_score = pct.1 - pct.2,
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

murthy <- readRDS("files/murthy_diff.RDS")
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

celegans <- fread("files/celegans_sc_ciliated_oxygen_2022.txt") %>%
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

trachea <- fread("files/trachea_all_markers.txt") %>%
  filter(
    cluster == "Ciliated Cells", pct.1 - pct.2 > 0,
    p_val_adj < 0.001
  ) %>%
  hgncConverter("gene_name") %>%
  mutate(
    mean_score = normalization0(pct.1 - pct.2),
    #mean_score = normalization0(log(pct.1 - pct.2)),
    #mean_score = pct.1 - pct.2,
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

cengen <- fread("files/difex_cengen.csv") %>%
  filter(human_symbol != "") %>%
  filter(
    cluster == "Sensory", pct.1 - pct.2 > 0,
    p_val_adj < 0.001
  ) %>%
  hgncConverter("human_symbol") %>%
  mutate(
    #mean_score = normalization0(pct.1 - pct.2),
    mean_score = normalization0(log(pct.1 - pct.2)),
    source = "cengen"
  ) %>%
  rename("Gene_name" = "human_symbol") %>%
  dplyr::select(Gene_name, mean_score, source) %>%
  group_by(Gene_name) %>%
  mutate(mean_score = max(mean_score)) %>%
  distinct()
cengen$mean_score[cengen$mean_score == 0] <- 0.00001



d <- density(single_cell_scores$single_cell_score[single_cell_scores$single_cell_score != 0]) # returns the density data
jpeg("xxx.jpeg")
plot(d) # plots the results
dev.off()

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

write.table(single_cell_scores, "single_cell_scores_16.10.23.txt", sep = "\t", row.names = FALSE, quote = FALSE)


## Comparative genomics ----

pwd <- "/home/kaplanlab/Desktop/Desktop_Items/MUSTAFA_PIR/BLAST_Cluster/Homo_sapiens/blast_folder2/results/"

newscores_50 <- filtering(pwd, 0.001, 50)

gower_dist_50 <- daisy(newscores_50[, 2:73], metric = "gower")
cluster_50 <- hclust(gower_dist_50, method = "ward.D2")

tree_50 <- cutree(cluster_50, k = 40)
aa_50 <- cbind(newscores_50[, c(1, 2)], tree_50) # Genes with clusters

aa_50_temp <- aa_50
colnames(aa_50_temp)[3] <- "cluster_num"

a <- as.data.frame(cbind(newscores_50, tree_50))
a[, 2:73] <- lapply(a[, 2:73], function(x) as.numeric(as.character(x)))


# dfx<-aa_50
# colnames(dfx)[3]<-"cluster_num"
# scores<-data.frame(cluster_num = c(31,37,8,2,17,39,16,3,30,10,23,24,1,4:7,9,11:15,18:22,25:29,32:36,38,40),
#                    score2 = c(1,0.9,0.6,0.6,0.5,0.5,0.4,0.3,0.2,0.1,0.1,0.1,rep(0,28)))
# **********************************************
dfx<-aa_50
colnames(dfx)[3]<-"cluster_num"
scores<-data.frame(cluster_num = c(31,37,1:30,32:36,38:40),
                   score2 = c(1,0.9,rep(0,38)))
# **********************************************
dfx<-hgncConverter(dfx, "Gene_name")
dfx<-left_join(dfx, scores, by = "cluster_num")
dfx$phylogenetic_scores<-normalization(dfx$score2)
cluster_scores<-left_join(hgnc_names[,2], dfx[,c(1,5)], by = c("Approved symbol" = "Gene_name"))
cluster_scores$phylogenetic_scores[is.na(cluster_scores$phylogenetic_scores)]<-0
cluster_scores<-cluster_scores %>%
  group_by(`Approved symbol`) %>%
  summarise_all(max)
colnames(cluster_scores)[1]<-"Gene_name"
write.table(cluster_scores, "cluster_scores.txt", sep = "\t", row.names = FALSE, quote = FALSE)


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

write.table(motif_scores2, "motif_scores.txt", sep = "\t", row.names = FALSE, quote = FALSE)



## Protein atlas scores ----

#source("proteinAtlas_webcrawling.R")

proatlas2 <- fread("files/protein_atlas_2022.txt")
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


## Statistial analysis ----

interaction_scores2 <- fread("protein_interaction_scores.txt")
cluster_scores2 <- fread("cluster_scores.txt")
motif_scores2 <- fread("motif_scores.txt")
proatlas_scores2 <- fread("proatlas_score.txt")
sc_scores2 <- fread("single_cell_scores_16.10.23.txt")
sc_scores2 <- fread("single_cell_scores.txt")


all_scores <- left_join(hgnc_names[, 2], interaction_scores2, by = c("Approved symbol" = "Gene_name"))
colnames(all_scores)[1] <- "Gene_name"
colnames(cluster_scores2)[1] <- "Gene_name"
all_scores <- left_join(all_scores, cluster_scores2, by = "Gene_name")
all_scores <- left_join(all_scores, motif_scores2, by = "Gene_name")
all_scores <- left_join(all_scores, proatlas_scores2, by = "Gene_name")
all_scores <- left_join(all_scores, sc_scores2, by = "Gene_name")

##****
ciliacarta<-read_xlsx("ciliacarta_2.xlsx")
ciliacarta<-ciliacarta[,c(2,7)]
ciliacarta<-unique(ciliacarta)
colnames(ciliacarta)<-c("Gene_name","ciliacarta")

ciliacarta<-ciliacarta[!duplicated(ciliacarta[,"Gene_name"]),]

ciliacarta<-hgncConverter(ciliacarta, "Gene_name")
ciliacarta<-unique(ciliacarta)

ciliacarta<-ciliacarta %>%
  group_by(Gene_name) %>%
  mutate(ciliacarta = mean(ciliacarta))
ciliacarta<-unique(ciliacarta)

all_scores <- left_join(all_scores, ciliacarta, by = "Gene_name")
all_scores1 <- rbind(
  semi_join(all_scores, ciliatest, by = "Gene_name"),
  semi_join(all_scores, negativeTest, by = "Gene_name")
)

##****
all_scores1<-all_scores1[!is.na(all_scores1$ciliacarta),]
all_scores1[is.na(all_scores1)] <- 0

# all_scores1 <- rbind(
#   semi_join(all_scores, ciliatest, by = "Gene_name"),
#   semi_join(all_scores, negativeTest, by = "Gene_name")
# )
all_scores1$total_score_sums <- rowSums(all_scores1[, 2:6])
all_scores1$is.ciliary <- ifelse(all_scores1$Gene_name %in% ciliatest$Gene_name, 1, 0)
all_scores1<-unique(all_scores1)

write.table(all_scores1, "score_table_roc_2.txt", row.names = FALSE, quote = FALSE, sep = "\t")

all_scores1<-fread("score_table_roc.txt")
#all_scores1$is.ciliary<-as.factor(all_scores1$is.ciliary)

### Determine thresholds

# th_interactions<-limiter(all_scores1, 2, 7) # 0.66
# th_phylogenetic<-limiter(all_scores1, 3, 7) # 0.61
# th_motif<-limiter(all_scores1, 4, 7) # 0.20
# th_pa<-limiter(all_scores1, 5, 7) # 0.25
# th_sc<-limiter(all_scores1, 6, 7) # 0.15
# 
# 
# for (i in 2:6) {
#   x <- classifier(all_scores1, i, "is.ciliary")
#   print(x)
# }

# Determine weight vector
df<-all_scores1

fn <- function(v) {
  pred <- df[[2]]*v[1]+df[[3]]*v[2]+df[[4]]*v[3]+df[[5]]*v[4]+df[[6]]*v[5]
  k <- df[[7]] - pred
  
  return <- sum(abs(k))
}

opt<-optim(c(0.5, 0.5, 0.5, 0.5, 0.5), fn)
opt<-optim(c(1, 1, 1, 1, 1), fn)

weight<-opt$par
all_scores1$total_score_sums <- rowSums(all_scores1[, 2:6])
ws <- as.data.frame(mapply(`*`, all_scores1[, 2:6], weight))
all_scores1$total_score_weighted <- rowSums(ws)

# *************************************** #
all_scores2<-all_scores
all_scores2$total_score_sums <- rowSums(all_scores2[, 2:6])
ws <- as.data.frame(mapply(`*`, all_scores2[, 2:6], weight))
all_scores2$total_score_weighted <- rowSums(ws)

all_scores2<-unique(all_scores2)

# for (i in c(2:6, 8:9)) {
#   x <- classifier2(all_scores1, i, "is.ciliary")
#   print(x)
# }

aaa<-limiter(all_scores1, 2, 7) # 0.58
aaa<-limiter(all_scores1, 3, 7) # 0.90
aaa<-limiter(all_scores1, 4, 7) # 0.20
aaa<-limiter(all_scores1, 5, 7) # 0.10
aaa<-limiter(all_scores1, 6, 7) # 0.08
aaa<-limiter(all_scores1, 8, 7) # 0.65
aaa<-limiter(all_scores1, 9, 7) # 0.41



## Figures ----

library(pROC)
library(purrr)
library(dplyr)
library(ggplot2)

## 1
rocs <- roc(is.ciliary ~ interaction_scores + phylogenetic_scores + motif_score +
              pa_score + single_cell_score + total_score_weighted, data = all_scores1)
#ggroc(rocs)

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
               color="darkgrey", linetype="dashed") + theme_classic()


## 2

rocs <- roc(is.ciliary ~ interaction_scores + phylogenetic_scores + motif_score +
              pa_score + single_cell_score + ciliacarta + total_score_sums, data = all_scores1)

lbls<-c("Interactions", "Comparative genomics", "TF binding", "Text mining", "scRNA-seq", "CiliaCarta", "CilioGenics")
#ggroc(rocs)

# extract auc
rocs %>% 
  map(~tibble(AUC = .x$auc)) %>% 
  bind_rows(.id = "Methods") -> data.auc

rocs %>% 
  map(~tibble(AUC = .x$auc)) %>% 
  rbindlist() -> data.auc

# generate labels labels
data.auc$Methods<-lbls
data.auc %>% 
  mutate(label_long=paste0(Methods," , AUC = ",paste(round(AUC,2))),
         label_AUC=paste0("AUC = ",paste(round(AUC,2)))) -> data.labels

# plot on a single plot with AUC in labels
ggroc(rocs, size = 4) +
  scale_color_discrete(labels=data.labels$label_long) +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1),
               color="darkgrey", linetype="dashed") +
  labs(color='Methods') + theme_bw(base_size = 30) +
  theme(legend.key.size = unit(2, "cm"))
  

ggsave("roc_all_bw_ccarta-yy.jpg", width = 21, height = 14, dpi = 300) 


### Upset ----

library(ComplexUpset)

interaction_scores2 <- fread("protein_interaction_scores.txt")
cluster_scores2 <- fread("cluster_scores.txt")
motif_scores2 <- fread("motif_scores.txt")
proatlas_scores2 <- fread("proatlas_score.txt")
sc_scores2 <- fread("single_cell_scores.txt")

# all<-interaction_scores2[,1]
# all<-lapply(dfs, function(x) left_join(all, get(x), by = "Gene_name"))
# 
# all<-cbind(interaction_scores2, sc_scores2, motif_scores2, cluster_scores2, proatlas_scores2)
# 
# colnames(interactions)[1]<-"Gene_name"

dfs<-c("interaction_scores2", "sc_scores2", "motif_scores2", "cluster_scores2", "proatlas_scores2")

for(df in dfs) {
  assign(df, setNames(get(df),  c("Gene_name","score")))
}

th_interactions<-limiter(all_scores1, 2, 7) # 0.57
th_phylogenetic<-limiter(all_scores1, 3, 7) # 0.89
th_motif<-limiter(all_scores1, 4, 7) # 0.20
th_pa<-limiter(all_scores1, 5, 7) # 0.25
th_sc<-limiter(all_scores1, 6, 7) # 0.08
th_total<-limiter(all_scores1, 8, 7) #0.65

all_methods<-left_join(interaction_scores2, sc_scores2, by = "Gene_name") %>%
  left_join(motif_scores2, by = "Gene_name") %>%
  left_join(cluster_scores2, by = "Gene_name") %>%
  left_join(proatlas_scores2, by = "Gene_name")

colnames(all_methods)<-c("Gene_name", "Interactions", "Single_cell", "TF_gene_interactions", "Comparative_genomics", 
                         "Protein_atlas")
all_methods<-all_methods %>%
  mutate(Interactions = ifelse(Interactions >= 0.66, 1, 0),
         Single_cell = ifelse(Single_cell >= 0.15, 1, 0),
         TF_gene_interactions = ifelse(TF_gene_interactions >= 0.20, 1, 0),
         Comparative_genomics = ifelse(Comparative_genomics >= 0.61, 1, 0),
         Protein_atlas = ifelse(Protein_atlas >= 0.25, 1, 0),
         ciliary = ifelse(Gene_name %in% ciliaryGenes$Gene_name, 1, 0),
         negative = ifelse(Gene_name %in% negativeCiliary$Gene_name, 1, 0)) %>%
  mutate(Type = ifelse(ciliary == 1, "Ciliary", ifelse(negative == 1, "Negative", "Candidate")))


all_methods<-all_scores2 %>%
  filter(total_score_sums >= 1.40)

colnames(all_methods)[c(1:6)]<-c("Gene_name", "Interactions", "Single_cell", "TF_gene_interactions", "Comparative_genomics", 
                         "Protein_atlas")

all_methods<-all_methods %>%
  mutate(Interactions = ifelse(Interactions >= 0.66, 1, 0),
         Single_cell = ifelse(Single_cell >= 0.15, 1, 0),
         TF_gene_interactions = ifelse(TF_gene_interactions >= 0.20, 1, 0),
         Comparative_genomics = ifelse(Comparative_genomics >= 0.61, 1, 0),
         Protein_atlas = ifelse(Protein_atlas >= 0.25, 1, 0),
         ciliary = ifelse(Gene_name %in% ciliaryGenes$Gene_name, 1, 0),
         negative = ifelse(Gene_name %in% negativeCiliary$Gene_name, 1, 0)) %>%
  mutate(Type = ifelse(ciliary == 1, "Ciliary", ifelse(negative == 1, "Negative", "Candidate")))

# all_methods<-rbind(interaction_scores2, sc_scores2, motif_scores2, cluster_scores2, proatlas_scores2) %>%
#   dplyr::select(Gene_name) %>%
#   distinct() %>% 
#   mutate(Interactions = ifelse(Gene_name %in% interaction_scores2$Gene_name, 1, 0),
#          Single_cell = ifelse(Gene_name %in% sc_scores2$Gene_name, 1, 0),
#          TF_gene_interactions = ifelse(Gene_name %in% motif_scores2$Gene_name, 1, 0),
#          Comparative_genomics = ifelse(Gene_name %in% cluster_scores2$Gene_name, 1, 0),
#          Protein_atlas = ifelse(Gene_name %in% proatlas_scores2$Gene_name, 1, 0),
#          ciliary = ifelse(Gene_name %in% ciliaryGenes$Gene_name, 1, 0),
#          negative = ifelse(Gene_name %in% negativeCiliary$Gene_name, 1, 0)) %>%
#   mutate(is.ciliary = ifelse(ciliary == 1, "Ciliary", ifelse(negative == 1, "Negative", "Unknown")))

all_methods2<-all_methods %>%
  filter_at(vars(-Gene_name, -ciliary, -negative, -Type), any_vars(. > 0)) %>%
  unique()

all_sc2<-all_sc %>% 
  filter(ciliary == 1 | negative == 1)

all_sc3<-all_sc %>%
  filter(negative == 0)

col<-brewer.pal(8, "Pastel2")

upset(all_methods, colnames(all_methods)[2:6],
      base_annotations=list(
        'Intersection'=list(
          #text_mapping=aes(label="hey"),
          counts=TRUE,
          aes=aes(x=intersection, fill=Type),
          geom=list(
            geom_bar(stat='count', position='stack', na.rm=TRUE),
            geom_text(
              aes(
                label=!!aes_percentage(relative_to='intersection')
              ),
              stat='count',
              size = 3, hjust = 0.5, vjust = 1, position = "stack",
              check_overlap = TRUE
              #position=position_fill(vjust = 1)
            ), scale_fill_manual(values = c("Ciliary" = col[2], "Candidate" = col[3], "Negative" = col[1]))
          )
        ) 
      ),
      set_sizes = upset_set_size(
        geom = geom_bar(aes(fill=Type)),
        position='right'
      ) + scale_fill_manual(values = c("Ciliary" = col[2], "Candidate" = col[3], "Negative" = col[1])),
      width_ratio=0.1, sort_intersections_by=c('degree','ratio'),
      guides='over'
)

ggsave("upset_sections_noOverlap.jpg", width = 18, height = 10, dpi = 300) 

library(ggrepel)
upset(all_methods, colnames(all_methods)[2:6],
      base_annotations=list(
        'Intersection'=list(
          #text_mapping=aes(label="hey"),
          counts=TRUE,
          aes=aes(x=intersection, fill=Type),
          geom=list(
            geom_bar(stat='count', position='stack', na.rm=TRUE),
            geom_text_repel(
              aes(
                label=!!aes_percentage(relative_to='intersection')
              ),
              stat='count',
              size = 3, hjust = 0.5, vjust = 1, position = "stack",
              max.overlaps = 45
              #position=position_fill(vjust = 1)
            ), scale_fill_manual(values = c("Ciliary" = col[2], "Candidate" = col[3], "Negative" = col[1]))
          )
        ) 
      ),
      set_sizes = upset_set_size(
        geom = geom_bar(aes(fill=Type)),
        position='right'
      ) + scale_fill_manual(values = c("Ciliary" = col[2], "Candidate" = col[3], "Negative" = col[1])),
      width_ratio=0.1, sort_intersections_by=c('degree','ratio'),
      guides='over'
)

ggsave("upset_sections_repelled.jpg", width = 18, height = 10, dpi = 300) 

upset(all_methods2, colnames(all_methods2)[2:6],
      base_annotations=list(
        'Intersection size'=intersection_size(
          counts=TRUE,
          mapping=aes(fill=is.ciliary),
          text_mapping=aes(label=paste0(round(
            !!get_size_mode('exclusive_intersection')/!!get_size_mode('inclusive_union') * 100
          ), '%')),
        )
      ),
      set_sizes = upset_set_size(
        geom = geom_bar(aes(fill=is.ciliary)),
        position='right'
      ),
      width_ratio=0.1)


### Interaction network ----

write.table(yy6, "table_network.txt", row.names = FALSE, quote = FALSE, sep = "\t")

library(GGally)

# Build network object

ciliarynet<-yy6[yy6$Interactor_A == "CC2D2A" | yy6$Interactor_B == "IFT74" | 
                  yy6$Interactor_A == "IFT74" | yy6$Interactor_B == "CC2D2A",]

glist<-c("WDR31", "WDR54", "WDR38", "ZC2HC1A", "ZNF474")
ciliarynet<-yy6[yy6$Interactor_A %in% glist | yy6$Interactor_B %in% glist,]

glist2<-ciliaryGenes$Gene_name
ciliarynet<-yy6[yy6$Interactor_A %in% glist2 | yy6$Interactor_B %in% glist2,]



net<-unique(c(ciliarynet$Interactor_A, ciliarynet$Interactor_B))

netnet<-yy6[yy6$Interactor_A %in% net & yy6$Interactor_B %in% net,]

ciliarynet<-unique(rbind(ciliarynet, netnet))

ciliarynet$is.ciliary<-ifelse(ciliarynet$Interactor_B %in% ciliaryGenes$Gene_name, "Ciliary", 
                              ifelse(ciliarynet$Interactor_B %in% negativeCiliary$Gene_name, "Non-ciliary", "Unknown"))

cilsimpnet<-unique(ciliarynet[,1:2])
cilsimpnet<-cilsimpnet[cilsimpnet$Interactor_A != cilsimpnet$Interactor_B,]

library(network)
nw <- network(cilsimpnet, directed = TRUE, matrix.type = "edgelist", loops = TRUE)

library(GGally)

candlist<-all_scores2$Gene_name[all_scores2$total_score_sums >= 1.2]

vertexnames<-network.vertex.names(nw)
nw %v% "is.ciliary" = ifelse(vertexnames %in% ciliaryGenes$Gene_name, "Ciliary", 
                             ifelse(vertexnames %in% negativeCiliary$Gene_name, "Non-ciliary",
                                    ifelse(vertexnames %in% candlist, "Candidate", "Unknown")))

set.seed(2)
ggnet2(nw, color = "is.ciliary", node.size = 5, label = c("CC2D2A", "IFT74"), label.size = 5,
       palette = c("Ciliary" = "steelblue", "Non-ciliary" = "tomato", "Unknown" = "grey70", "Candidate" = "yellow2"),
       mode = "kamadakawai") +
  guides(color=guide_legend(title="Type", override.aes = list(size=7))) +
  theme(legend.title = element_text(size=15),
        legend.text = element_text(size=12))

ggsave("interaction_plot_IFT74_CC2D2A_2-2.jpg", width = 12, height = 8, dpi = 300) 

nw %v% "size" = ifelse(vertexnames %in% glist, 5, 3)

set.seed(29)
ggnet2(nw, color = "is.ciliary", node.size = 5, label = TRUE, label.size = "size",
       palette = c("Ciliary" = "steelblue", "Non-ciliary" = "tomato", "Unknown" = "grey70", "Candidate" = "yellow2"),
       mode = "fruchtermanreingold") +
  guides(color=guide_legend(title="Type", override.aes = list(size=6))) +
  theme(legend.title = element_text(size=15),
        legend.text = element_text(size=12))


ggsave("interaction_plot_candidates_all2-2.jpg", width = 12, height = 8, dpi = 300) 

ggnet2(nw, color = "is.ciliary", node.size = 5, label = c("ARL13B", "IFT74"), label.size = 4,
       palette = c("Ciliary" = col[3], "Non-ciliary" = col[1], "Unknown" = col[2]),
       mode = "kamadakawai")

library(RColorBrewer)
col<-brewer.pal(8, "Set3")[c(1,4,9)]




### Interaction tables ----

library(formattable)
library(htmlwidgets)
library(webshot2)
library(stringr)
library(tidyr)

pub<-fread("publications.txt")
pub<-pub[!is.na(pub$Gene_name),]

pub<-hgncConverter(pub, "Gene_name")
pub2<-pub %>%
  group_by(Gene_name) %>%
  count(Gene_name)
ciliogenics<-all_scores2
colnames(ciliogenics)[7]<-"total_score"
ciliogenics1<-ciliogenics[order(ciliogenics$total_score, decreasing = TRUE),]
#ciliogenics1$publication<-ifelse(ciliogenics1$Gene_name %in% pub$Gene_name, 1, 0)
ciliogenics1<-left_join(ciliogenics1, pub2, by = "Gene_name")
colnames(ciliogenics1)[9]<-"Publication"

#yy7<-yy6 %>% left_join(ciliogenics1[,c(1,9)], by = "Gene_name")

int_table<-yy6[1:20,]

int_tbl<-interaction_scores2
colnames(int_tbl)[2]<-"score"
int_tbl<-int_tbl[order(int_tbl$score, decreasing = TRUE),]
colnames(pub2)[2]<-"Publications"

int_tbl<-int_tbl %>% 
  left_join(pub2, by = "Gene_name")

int_tbl2<-int_tbl %>% 
  left_join(pub2, by = "Gene_name") %>%
  filter(!(Gene_name %in% ciliaryGenes$Gene_name)) %>%
  arrange(desc(score)) %>%
  slice(1:20) %>%
  left_join(yy6[,1:2], by = c("Gene_name" = "Interactor_A")) %>%
  left_join(yy6[,1:2], by = c("Gene_name" = "Interactor_B")) %>%
  pivot_longer(cols = Interactor_B:Interactor_A, values_to = 'Interactor_B') %>%
  dplyr::select(-name) %>%
  distinct() %>%
  group_by(Gene_name, score) %>%
  slice(1:5) %>%
  summarize(Interactors = str_c(Interactor_B, collapse = ", ")) %>%
  mutate(`scRNA-seq` = ifelse(Gene_name %in% sc_scores2$Gene_name[sc_scores2$single_cell_score >= 0.08], "Yes", "No"),
         `Comparative Genomics` = ifelse(Gene_name %in% cluster_scores2$`Approved symbol`[cluster_scores2$phylogenetic_scores >= 0.61], "Yes", "No"),
         `Publications` = int_tbl$Publications[match(Gene_name, int_tbl$Gene_name)]) %>%
  rename("Gene Name" = "Gene_name",
         "Interaction score" = "score",
         "Number of publications" = "Publications")

proat_loc<-progenes[progenes$Cilia_localization == 1,]
int_tbl_loc<-int_tbl2
int_tbl_loc$`Cilia localization`<-ifelse(int_tbl_loc$`Gene Name` %in% proat_loc$Gene_name, "Yes", "NA")
ddd<-c("Basal body","Centrosome","Centrosome","Centrosome","Centrosome","Centrosome","Basal body","Centrosome","Cilia base","NA",
       "NA","NA","NA","Centrosome","NA","NA","Cilia","Centrosome","NA","NA")
int_tbl_loc$`Cilia localization`<-ddd


customGreen0 = "#DeF7E9"

customGreen = "#71CA97"

customRed = "#ff7f7f"

tbl<-formattable(int_tbl2,
                 align =c("l","c","l","c","c","c"),
                 list(`Gene Name` = formatter(
                   "span", style = ~ style(color = "grey",font.weight = "bold")),
                   `Interaction score`= color_tile(customGreen0, customGreen),
                   `scRNA-seq` = formatter("span",
                                          style = x ~ style(color = ifelse(x == "Yes", "green", "red")),
                                          x ~ icontext(ifelse(x == "Yes", "ok", "remove"), ifelse(x == "Yes", "Yes", "No"))),
                   `Comparative Genomics` = formatter("span",
                                           style = x ~ style(color = ifelse(x == "Yes", "green", "red")),
                                           x ~ icontext(ifelse(x == "Yes", "ok", "remove"), ifelse(x == "Yes", "Yes", "No"))),
                   `Number of publications` = color_tile(customGreen0, customGreen)
                   ), format = "f", digits = 3)
tbl
                   
w <- as.htmlwidget(tbl)
saveWidget(w, "table.html", selfcontained = TRUE)

webshot2::webshot(url = "table.html", file = "table.png", 
                  vwidth = 1500, vheight = 2000)



wrap_in <- function(.df, column, word, tag){
  class<-""
  if(grepl("\\.", tag)) {
    class <- sub(".+?\\.(.+)", " class='\\1'", tag)
    tag <- sub("\\..+", "", tag)
  }
  .df[[column]] <-  gsub(sprintf("\\b(%s)\\b", paste0(word,collapse="|")), sprintf("<%1$s%2$s>\\1</%1$s>", tag, class), .df[[column]])
  .df
}


int_tbl_loc %>% wrap_in( "Interactors", ciliaryGenes$Gene_name, "span.highlight") %>% 
  formattable(align =c("l","c","l","c","c","c"),
              list(`Gene Name` = formatter(
                "span", style = ~ style(color = "grey",font.weight = "bold")),
                `Interaction score`= color_tile(customGreen0, customGreen),
                `scRNA-seq` = formatter("span",
                                        style = x ~ style(color = ifelse(x == "Yes", "green", "red")),
                                        x ~ icontext(ifelse(x == "Yes", "ok", "remove"), ifelse(x == "Yes", "Yes", "No"))),
                `Comparative Genomics` = formatter("span",
                                                   style = x ~ style(color = ifelse(x == "Yes", "green", "red")),
                                                   x ~ icontext(ifelse(x == "Yes", "ok", "remove"), ifelse(x == "Yes", "Yes", "No"))),
                `Number of publications` = color_tile(customGreen0, customGreen),
                `Cilia localization` = formatter("span",
                                        style = x ~ style(color = ifelse(x != "NA", "green", "red")),
                                        x ~ icontext(ifelse(x != "NA", "ok", "remove"), ifelse(x != "NA", x, "NA")))
              ), format = "f", digits = 3) -> tbl

w <- as.htmlwidget(tbl)
w1<-w %>% htmlwidgets::prependContent(htmltools::tags$style("span.highlight {color:red;font-weight: 700;}"))

saveWidget(w1, "table_2_loc.html", selfcontained = TRUE)

webshot2::webshot(url = "tablex.html", file = "table.pdf", delay = 2)



### Comparative genomics tables ----

library(gprofiler2)

lgost<-list()
for(i in 1:length(unique(aa_50$tree_50))){
  gost<-gost(query = aa_50$Gene_name[aa_50$tree_50 == i])
  names(gost)<-paste0("cluster_", i)
  lgost<-c(lgost, gost[1])
}

gotable<-function(cluster_number) {
  a<-lgost[[cluster_number]]
  a<-rbind(head(a[a$source == "GO:BP",]),
           head(a[a$source == "GO:MF",]),
           head(a[a$source == "GO:CC",]))
  
  a$query<-paste0("cluster_", cluster_number)
  golist<-a$term_id[1:5]
  
  gost<-gost(query = aa_50$Gene_name[aa_50$tree_50 == cluster_number])
  tt<-publish_gosttable(gost, highlight_terms = golist)
  
  ggsave(paste0("gosttable_cluster", cluster_number, ".jpg"), width = 12, height = 8, dpi = 300) 
  
}

gotable(31)
gotable(37)
gotable(1)
gotable(24)


### Heatmap- Comparative genomics ----

library(iheatmapr)

nscores2<-newscores_50
nscores2[duplicated(nscores2$Gene_name),1]<-paste0(nscores2[duplicated(nscores2$Gene_name),1], "_2")
nscores2<-nscores2[!is.na(nscores2$Gene_name),]
row.names(nscores2)<-nscores2$Gene_name

aa_50

aa<-aa_50
aa<-aa[,-2]
colnames(aa)[2]<-"cluster_number"
#aa<-hgncConverter(aa, "Gene_name")

species<-fread("species_app.txt")
species<-species[,-1]

species<-species[which(species$X1 %in% colnames(newscores_50[,2:73])),]
anot<-data.frame(Class = species[,2], Organisms = rep(c("Ciliary", "Nonciliary"), c(44,28)))
row.names(anot)<-species$X1
colnames(anot)[1]<-"Class"

my_colour = list(Organisms = c(Ciliary = "firebrick3", Nonciliary = "dodgerblue3"),
                 Class = c(Animals = "firebrick3", Fungi = "dodgerblue3", Protists = "darkgrey", Plants = "chartreuse", Other = "ghostwhite", Bacteria = "gray0"))

cnumber<-31
a<-as.matrix(nscores2[which(nscores2$cluster_number == cnumber),2:73])
rownames(a)<-nscores2$Gene_name[which(nscores2$cluster_number == cnumber)]
a

main_heatmap(a, layout = list(paper_bgcolor='transparent'),
             tooltip = setup_tooltip_options(prepend_row = "Gene: ", prepend_col = "Organism: "))%>%
  add_col_labels(size = 0.46,font = list(family = c("open_sansregular"), size = 14), textangle=90,
                 tickvals = c(1:length(colnames(a))))%>%
  add_col_annotation(annotation=anot, side="top", size = 0.1) %>%
  modify_layout(list(margin = list(l = 80)))


### Motif network ----

ciliarynet<-yy6[yy6$Interactor_A == "CC2D2A" | yy6$Interactor_B == "IFT74" | 
                  yy6$Interactor_A == "IFT74" | yy6$Interactor_B == "CC2D2A",]

glist<-c("WDR31", "WDR54", "WDR38", "ZC2HC1A", "ZNF474","BASP1")
glist<-c("IFT88", "TMEM231", "NEK10",  "WDR38", "WDR54", "ZNF474", "WDR38", "BASP1", "ZC2HC1A", "TMEM190", "TMEM145", "WDR31")
motifnet<-pan_motif[pan_motif$Gene_name %in% glist,]

glist2<-ciliaryGenes$Gene_name
motifnet2<-pan_motif[pan_motif$Gene_name %in% glist2,]



net<-unique(c(pan_motif$Gene_name))

netnet<-pan_motif[pan_motif$Gene_name %in% net,]



motifnet$is.ciliary<-ifelse(motifnet$Gene_name %in% ciliaryGenes$Gene_name, "Ciliary", 
                              ifelse(motifnet$Gene_name %in% negativeCiliary$Gene_name, "Non-ciliary", "Unknown"))

motsimpnet<-unique(motifnet[,c(1,3)])


library(network)
nw <- network(motsimpnet, directed = TRUE, matrix.type = "edgelist", loops = TRUE)

library(GGally)
ggnet2(nw)

candlist<-c("WDR31","ZC2HC1A","BASP1")
cilmotiflist<-c("RFX2", "RFX3", "MYB", "GLIS3", "JAZF1", "FOXJ1", "SOX5", "TOX")

vertexnames<-network.vertex.names(nw)
nw %v% "is.ciliary" = ifelse(vertexnames %in% ciliaryGenes$Gene_name | vertexnames %in% cilmotiflist, "Ciliary", 
                             ifelse(vertexnames %in% negativeCiliary$Gene_name, "Non-ciliary",
                                    ifelse(vertexnames %in% candlist, "Candidate", "Non-ciliary")))

nw %v% "size" = ifelse(vertexnames %in% glist, 7, 5)

set.seed(2)
ggnet2(nw, color = "is.ciliary", node.size = 8, label = TRUE, label.size = "size",
       palette = c("Ciliary" = "steelblue", "Non-ciliary" = "tomato", "Candidate" = "yellow2"),
       mode = "kamadakawai") +
  guides(color=guide_legend(title="Type", override.aes = list(size=8))) +
  theme(legend.title = element_text(size=15),
        legend.text = element_text(size=12))


ggsave("motif_plot.jpg", width = 12, height = 10, dpi = 300) 


### Motif tables ----
motif_scores3<-motif_scores2
colnames(motif_scores3)[2]<-"score"

mot_table<-left_join(pan_motif, motif_scores3, by = "Gene_name") %>%
  dplyr::select(-2) %>%
  filter(!(Gene_name %in% ciliaryGenes$Gene_name)) %>%
  group_by(Gene_name, score) %>%
  slice(1:5) %>%
  summarize(TFs = str_c(motif, collapse = ", ")) %>%
  na.omit() %>%
  ungroup() %>%
  arrange(desc(score)) %>%
  slice(1:10) %>%
  mutate(`scRNA-seq` = ifelse(Gene_name %in% sc_scores2$Gene_name[sc_scores2$single_cell_score >= 0.08], "Yes", "No"),
         `Comparative Genomics` = ifelse(Gene_name %in% cluster_scores2$`Approved symbol`[cluster_scores2$phylogenetic_scores >= 0.61], "Yes", "No"),
         `Publications` = int_tbl$Publications[match(Gene_name, int_tbl$Gene_name)]) %>%
  rename("Gene Name" = "Gene_name")

mot_table1<-left_join(pan_motif, motif_scores3, by = "Gene_name") %>%
  dplyr::select(-2) %>%
  filter(!(Gene_name %in% ciliaryGenes$Gene_name)) %>%
  group_by(Gene_name, score) %>%
  slice(1:5) %>%
  summarize(TFs = str_c(motif, collapse = ", ")) %>%
  na.omit() %>%
  ungroup() %>%
  arrange(score) %>%
  slice(1:5) %>%
  mutate(`scRNA-seq` = ifelse(Gene_name %in% sc_scores2$Gene_name[sc_scores2$single_cell_score >= 0.08], "Yes", "No"),
         `Comparative Genomics` = ifelse(Gene_name %in% cluster_scores2$`Approved symbol`[cluster_scores2$phylogenetic_scores >= 0.61], "Yes", "No"),
         `Publications` = int_tbl$Publications[match(Gene_name, int_tbl$Gene_name)]) %>%
  rename("Gene Name" = "Gene_name")

genelist<-c("IFT88", "TMEM231", "NEK10",  "WDR38", "WDR54", "ZNF474", "WDR38", "BASP1", "ZC2HC1A", "TMEM190", "TMEM145", "WDR31")
mot_table2<-left_join(pan_motif, motif_scores3, by = "Gene_name") %>%
  dplyr::select(-2) %>%
  filter(!(Gene_name %in% ciliaryGenes$Gene_name)) %>%
  group_by(Gene_name, score) %>%
  slice(1:5) %>%
  summarize(TFs = str_c(motif, collapse = ", ")) %>%
  na.omit() %>%
  ungroup() %>%
  filter(Gene_name %in% genelist) %>%
  mutate(`scRNA-seq` = ifelse(Gene_name %in% sc_scores2$Gene_name[sc_scores2$single_cell_score >= 0.08], "Yes", "No"),
         `Comparative Genomics` = ifelse(Gene_name %in% cluster_scores2$`Approved symbol`[cluster_scores2$phylogenetic_scores >= 0.61], "Yes", "No"),
         `Publications` = int_tbl$Publications[match(Gene_name, int_tbl$Gene_name)]) %>%
  rename("Gene Name" = "Gene_name")

mot_table<-rbind(mot_table, mot_table2, mot_table1)
mot_table$Publications[is.na(mot_table$Publications)]<-0

colnames(mot_table)[c(2,6)]<-c("TF binding score", "Number of publications")

proat_loc<-progenes[progenes$Cilia_localization == 1,]
mot_table_loc<-mot_table
mot_table_loc$`Cilia localization`<-ifelse(mot_table_loc$`Gene Name` %in% proat_loc$Gene_name, "Yes", "No")
ddd1<-c()

customGreen0 = "#DeF7E9"
customGreen = "#71CA97"
customRed = "#ff7f7f"

motif_list<-c("RFX2", "RFX3", "MYB", "GLIS3", "JAZF1", "FOXJ1", "SOX5", "TOX")

wrap_in <- function(.df, column, word, tag){
  class<-""
  if(grepl("\\.", tag)) {
    class <- sub(".+?\\.(.+)", " class='\\1'", tag)
    tag <- sub("\\..+", "", tag)
  }
  .df[[column]] <-  gsub(sprintf("\\b(%s)\\b", paste0(word,collapse="|")), sprintf("<%1$s%2$s>\\1</%1$s>", tag, class), .df[[column]])
  .df
}

mot_table_loc %>% wrap_in("TFs", motif_list, "span.highlight") %>% 
  formattable(align =c("l","c","l","c","c","c","c"),
              list(`Gene Name` = formatter(
                "span", style = ~ style(color = "grey",font.weight = "bold")),
                `TF binding score`= color_tile(customGreen0, customGreen),
                `scRNA-seq` = formatter("span",
                                        style = x ~ style(color = ifelse(x == "Yes", "green", "red")),
                                        x ~ icontext(ifelse(x == "Yes", "ok", "remove"), ifelse(x == "Yes", "Yes", "No"))),
                `Comparative Genomics` = formatter("span",
                                                   style = x ~ style(color = ifelse(x == "Yes", "green", "red")),
                                                   x ~ icontext(ifelse(x == "Yes", "ok", "remove"), ifelse(x == "Yes", "Yes", "No"))),
                `Number of publications` = color_tile(customGreen0, customGreen),
                `Cilia localization` = formatter("span",
                                        style = x ~ style(color = ifelse(x == "Yes", "green", "red")),
                                        x ~ icontext(ifelse(x == "Yes", "ok", "remove"), ifelse(x == "Yes", "Yes", "No")))
              ), format = "f", digits = 3) -> tbl

w <- as.htmlwidget(tbl)
w1<-w %>% htmlwidgets::prependContent(htmltools::tags$style("span.highlight {color:red;font-weight: 700;}"))

saveWidget(w1, "table_motif3.html", selfcontained = TRUE)

webshot2::webshot(url = "table_motif3.html", file = "table_motif3.pdf", delay = 2)



# Venn diagram

library(patchwork)
library(VennDiagram)
library(RColorBrewer)

mycol2 <- brewer.pal(3, "Pastel2")

graphics.off()
venn<-venn.diagram(
  x = list("CilioGenics" = ciliogenics$Gene_name[ciliogenics$phylogenetic_scores > 0.8],
           "Nevers et al." = yannis$Gene_name,
           "CLIME" = clime$Gene_name),
  lwd = 2,
  lty = 'blank',
  fill = mycol2,
  cex = 2,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  #cat.pos = c(-27, 27, 135),
  #cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1,
  filename = NULL
)

grid::grid.draw(venn)
ggsave(plot = venn, "venn_comparative_genomics.jpg", width = 12, height = 10, dpi = 300)


### Heatmap - compare all methods ----

alll<-all_scores2 %>% 
  dplyr::select(-total_score_weighted) %>%
  mutate(is.ciliary = ifelse(Gene_name %in% negativeCiliary$Gene_name, "Negative", 
                             ifelse(Gene_name %in% ciliaryGenes$Gene_name, "Ciliary", "Unknown"))) %>%
  group_by(Gene_name) %>%
  unique() %>%
  tibble::column_to_rownames("Gene_name")


hgnc_names1<-unique(hgnc_names[,2])

alll<-left_join(hgnc_names1, reyfman[,1:2], by = c("Approved symbol" = "Gene_name"))

alll<-left_join(alll, carraro[,1:2], by = c("Approved symbol" = "Gene_name"))
alll<-left_join(alll, habermann[,1:2], by = c("Approved symbol" = "Gene_name"))
alll<-left_join(alll, celegans[,1:2], by = c("Approved symbol" = "Gene_name"))
colnames(alll)<-c("Gene_name", "reyfman", "carraro", "habermann", "celegans")

interaction_scores<-interaction_scores2 %>%
  group_by(Gene_name) %>%
  mutate(interaction_scores = max(interaction_scores)) %>%
  unique()

cluster_scores<-cluster_scores2 %>%
  group_by(Gene_name) %>%
  mutate(phylogenetic_scores = max(phylogenetic_scores)) %>%
  unique()

motif_scores<-motif_scores2 %>%
  group_by(Gene_name) %>%
  mutate(motif_score = max(motif_score)) %>%
  unique()

proatlas_scores<-proatlas_scores2 %>%
  group_by(Gene_name) %>%
  mutate(pa_score = max(pa_score)) %>%
  unique()


alll<-left_join(alll, interaction_scores, by = "Gene_name")
alll<-left_join(alll, cluster_scores, by = "Gene_name")
alll<-left_join(alll, motif_scores, by = "Gene_name")
alll<-left_join(alll, proatlas_scores, by = "Gene_name")
alll<-left_join(alll, single_cell_scores, by = "Gene_name")

alll$is.ciliary<-ifelse(alll$Gene_name %in% negativeCiliary$Gene_name, "Negative", 
                        ifelse(alll$Gene_name %in% ciliaryGenes$Gene_name, "Ciliary", "Unknown"))

alll[is.na(alll),]<-0


alll<-alll %>% 
  group_by(Gene_name) %>%
  unique() %>%
  tibble::column_to_rownames("Gene_name")

alll$is.ciliary<-as.factor(alll$is.ciliary)

alll$total<-rowSums(alll[,6:9])

alll2<-alll %>%
  filter(is.ciliary %in% c("ciliary", "negative") | total > 1)

alll2<-rbind(alll2[alll2$is.ciliary == "ciliary",], alll2[alll2$is.ciliary == "negative",], alll2[alll2$is.ciliary == "unknown",])

anot<-data.frame(Type = alll2$is.ciliary, stringsAsFactors = TRUE)
rownames(anot)<-alll2$Gene_name
#rownames(cg1)<-cg1[,1]


an_col<-list(Type = c(ciliary = "royalblue1", negative = "steelblue1", unknown = "lightskyblue1"))
col<-c("lightskyblue", "royalblue4")
col<-brewer.pal(9, "GnBu")
col<-brewer.pal(9, "YlGn")
col<-brewer.pal(3, "Pastel2")

rownames(alll2)<-alll2$Gene_name
alll2<-as.data.frame(alll2)

pa<-pheatmap(alll2[,5:9], cluster_cols = FALSE, cluster_rows = FALSE, 
         show_rownames = FALSE, angle_col = 315, annotation_row = anot,
         annotation_colors = an_col, scale = "column")

ggsave("heatmap_all_methods.jpg", width = 9, height = 12, dpi = 300) 

d<-density(alll2$single_cell_score[alll2$single_cell_score != 0])
plot(d)

alll3<-alll2 %>% 
  mutate(interaction_scores = ifelse(interaction_scores > 0.4, 1, 0)) %>%
  mutate(single_cell_score = ifelse(single_cell_score > 0.2, 1, 0)) %>%
  mutate(phylogenetic_scores = ifelse(phylogenetic_scores > 0.8, 1, 0)) %>%
  mutate(motif_score = ifelse(motif_score > 0.5, 1, 0)) %>%
  mutate(pa_score = ifelse(pa_score > 0.8, 1, 0))

pheatmap(alll3[,6:10], cluster_cols = FALSE, cluster_rows = FALSE, 
         show_rownames = FALSE, angle_col = 315, annotation_row = anot,
         annotation_colors = an_col, color = col[c(2,5)])







alll$is.ciliary<-as.factor(alll$is.ciliary)

alll2<-alll %>%
  filter(is.ciliary %in% c("Ciliary", "Negative") | total_score_sums > 0.68)

alll2<-rbind(alll2[alll2$is.ciliary == "Ciliary",], alll2[alll2$is.ciliary == "Negative",], alll2[alll2$is.ciliary == "Unknown",])

anot<-data.frame(Type = alll2$is.ciliary, stringsAsFactors = TRUE)
rownames(anot)<-rownames(alll2)
#rownames(cg1)<-cg1[,1]

an_col<-list(Type = c(Ciliary = "royalblue1", Negative = "steelblue1", Unknown = "lightskyblue1"))
col<-c("lightskyblue", "royalblue4")
col<-brewer.pal(9, "GnBu")
col<-brewer.pal(9, "YlGn")
col<-brewer.pal(3, "Pastel2")

alll2<-as.data.frame(alll2)
alll2[alll2$is.ciliary == "Unknown",]<-alll2[order(alll2$total_score_sums[alll2$is.ciliary == "Unknown"], decreasing = TRUE),]
pa<-pheatmap(alll2[,1:5], cluster_cols = FALSE, cluster_rows = FALSE, 
         show_rownames = FALSE, angle_col = 315, annotation_row = anot,
         annotation_colors = an_col, scale = "column", fontsize = 9,
         labels_col = c("Interactions", "Comparative genomics", "TF binding", "Text mining", "scRNA-seq"),
         width = 10)

ggsave("heatmap_all_methods3.jpg", plot = pa, width = 5, height = 5, dpi = 300)

d<-density(alll2$single_cell_score[alll2$single_cell_score != 0])
plot(d)

alll3<-alll2 %>% 
  mutate(interaction_scores = ifelse(interaction_scores > 0.4, 1, 0)) %>%
  mutate(single_cell_score = ifelse(single_cell_score > 0.2, 1, 0)) %>%
  mutate(phylogenetic_scores = ifelse(phylogenetic_scores > 0.8, 1, 0)) %>%
  mutate(motif_score = ifelse(motif_score > 0.5, 1, 0)) %>%
  mutate(pa_score = ifelse(pa_score > 0.8, 1, 0))

pheatmap(alll3[,1:5], cluster_cols = FALSE, cluster_rows = FALSE, 
         show_rownames = FALSE, angle_col = 315, annotation_row = anot,
         annotation_colors = an_col)


### Protein atlas figure ----


proatlas_scores2


proatlas2 <- fread("files/protein_atlas_2022.txt")
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

# ************************************************
progenes$Centrosome <- ifelse(progenes$Centrosome == "YES", 1, 0)
progenes$Cilia_expression <- ifelse(progenes$Cilia_expression == "YES", 1, 0)
progenes$Flagella <- ifelse(progenes$Flagella == "YES", 1, 0)
progenes$Cilia_localization <- ifelse(progenes$Cilia_localization == "YES", 1, 0)
# ************************************************

progenes$Type<-ifelse(progenes$Gene_name %in% ciliaryGenes$Gene_name, "Ciliary", 
                      ifelse(progenes$Gene_name %in% negativeCiliary$Gene_name, "Negative", "Unknown"))

progenes<-progenes %>% full_join(allgenes, by = "Gene_name")
progenes$is.ciliary[is.na(progenes$is.ciliary)]<-2

progenes$Type[is.na(progenes$Type)]<-ifelse(progenes$is.ciliary[is.na(progenes$Type)] == 1, "Ciliary",
                                            ifelse(progenes$is.ciliary[is.na(progenes$Type)] == 0, "Negative", "Unknown"))
progenes[is.na(progenes)]<-0
progenes<-rbind(progenes[progenes$Type == "Ciliary",],
                progenes[progenes$Type == "Negative",],
                progenes[progenes$Type == "Unknown",])

plos<-read_xlsx("plos1_genes.xlsx") %>%
  hgncConverter("Gene_name")

plos$`Ivliev et al.`<-1

progenes<-left_join(progenes, plos, by = "Gene_name")
progenes[is.na(progenes)]<-0
progenes<-progenes[,-7]

anot<-data.frame(Type = progenes$Type, stringsAsFactors = TRUE)
rownames(anot)<-progenes$Gene_name
progenes<-as.data.frame(progenes)
rownames(progenes)<-progenes$Gene_name


an_col<-list(Type = c(Ciliary = "royalblue4", Negative = "royalblue1", Unknown = "skyblue1"))
an_col<-list(Type = c(Ciliary = "royalblue1", Negative = "steelblue1", Unknown = "lightskyblue1"))
col<-c("lightskyblue", "royalblue4")
col<-brewer.pal(9, "GnBu")
col<-brewer.pal(9, "YlGn")
col<-brewer.pal(8, "Pastel2")

ph<-pheatmap(progenes[,c(2:5,7)], cluster_cols = FALSE, cluster_rows = FALSE, 
             show_rownames = FALSE, angle_col = 315, annotation_row = anot,
             annotation_colors = an_col, color = col[c(2,4)])

ggsave("pa_comparison2.jpg", plot = ph, width = 4.5, height = 6.5, dpi = 300)


### SC comparison heatmap ----

yannis<-read_xlsx(paste0(pwd, "/files/SupplementaryTable2_edited.xlsx")) %>%
  rename("Gene_name" = "Gene Name") %>%
  mutate(yannis = 1)
clime<-fread(paste0(pwd,"/files/human_cilia_clime.txt"), select = "Gene Symbol") %>%
  rename("Gene_name" = "Gene Symbol") %>%
  mutate(clime = 1)


setwd("/home/kaplanlab/Desktop/Desktop_Items/MUSTAFA_PIR/ciliogenics")
cluster_scores<-fread("cluster_scores_2022_v3.txt")

comp<-left_join(cluster_scores2, clime, by = c("Approved symbol" = "Gene_name")) %>%
  left_join(yannis, by = c("Approved symbol" = "Gene_name")) %>%
  rename("Gene_name" = "Approved symbol")
comp[is.na(comp),]<-0

comp<-comp %>%
  filter_at(vars(-Gene_name), any_vars(. > 0)) %>%
  unique() %>%
  mutate(is.ciliary = ifelse(Gene_name %in% ciliaryGenes$Gene_name, "Ciliary", 
                             ifelse(Gene_name %in% negativeCiliary, "Negative", "Unknown")))


mycol2 <- brewer.pal(3, "Pastel2")

library(pheatmap)

allgenes<-unique(rbind(ciliaryGenes, negativeCiliary))
allgenes<-allgenes[!duplicated(allgenes$Gene_name),]
allgenes<-as.data.table(allgenes)
cg3<-left_join(allgenes[,1], comp, by = "Gene_name")
cg3<-cg3 %>% mutate(is.ciliary = ifelse(Gene_name %in% ciliaryGenes$Gene_name, "Ciliary", 
                                        ifelse(Gene_name %in% negativeCiliary$Gene_name, "Negative", "Unknown")))
cg3[is.na(cg3)]<-0

nocomp<-anti_join(comp, allgenes, by = "Gene_name")
cg3<-rbind(cg3, nocomp)

colnames(cg3)[2:5]<-c("CilioGenics", "CLIME", "Nevers et al.", "Type")

anot<-data.frame(Type = cg3$Type, stringsAsFactors = TRUE)
rownames(anot)<-cg3$Gene_name
cg3<-as.data.frame(cg3)
rownames(cg3)<-cg3$Gene_name


an_col<-list(Type = c(Ciliary = "royalblue4", Negative = "royalblue1", Unknown = "skyblue1"))
an_col<-list(Type = c(Ciliary = "royalblue1", Negative = "steelblue1", Unknown = "lightskyblue1"))
col<-c("lightskyblue", "royalblue4")
col<-brewer.pal(9, "GnBu")
col<-brewer.pal(9, "YlGn")
col<-brewer.pal(8, "Pastel2")

ph<-pheatmap(cg3[,2:4], cluster_cols = FALSE, cluster_rows = FALSE, 
             show_rownames = FALSE, angle_col = 315, annotation_row = anot,
             annotation_colors = an_col, color = col[c(2,4)])

ggsave("cg_comparison.jpg", plot = ph, width = 3, height = 6, dpi = 300) 


### SC comparison table ----

overlap<-calculate.overlap(list(`Carraro et al.` = carraro$Gene_name,
                                `Habermann et al.` = habermann$Gene_name,
                                `Reyfman et al.` = reyfman$Gene_name,
                                `Cao et al.` = celegans$Gene_name,
                                `Murthy et al.` = murthy$Gene_name))
overlap_list<-overlap$a31
sc_overlaps<-all_scores2 %>%
  filter(Gene_name %in% overlap_list) %>%
  mutate(Type = ifelse(Gene_name %in% ciliaryGenes$Gene_name, "Ciliary",
                       ifelse(Gene_name %in% negativeCiliary$Gene_name, "Negative", "Candidate")),
         phylogenetic_scores = ifelse(phylogenetic_scores > 0.8, "Yes", "No")) %>%
  rename(`Comparative Genomics` = "phylogenetic_scores")

sc_overlaps[,c(1,3,9)] %>% 
  formattable(align =c("l","c","l","c","c","c"),
              list(`Gene Name` = formatter(
                "span", style = ~ style(color = "grey",font.weight = "bold")),
                `Comparative Genomics` = formatter("span",
                                                   style = x ~ style(color = ifelse(x == "Yes", "green", "red")),
                                                   x ~ icontext(ifelse(x == "Yes", "ok", "remove"), ifelse(x == "Yes", "Yes", "No"))),
                `Type` = formatter("span",
                                   style = x ~ style(color = ifelse(x == "Ciliary", "green", ifelse(x == "Negative", "red", "blue"))),
                                   x ~ icontext(ifelse(x == "Ciliary", "ok", ifelse(x == "Negative", "remove", "thumbs-up")), 
                                                ifelse(x == "Ciliary", "Ciliary", ifelse(x == "Negative", "Negative", "Candidate"))))
              ), format = "f", digits = 3) -> tbl

w <- as.htmlwidget(tbl)
w1<-w %>% htmlwidgets::prependContent(htmltools::tags$style("span.highlight {color:red;font-weight: 700;}"))

saveWidget(w1, "table_sc_overlaps.html", selfcontained = TRUE)


## SC - Comparative genomics comparison table ----

genelist<-c("IFT88", "TMEM231", "NEK10", "WDR54", "ZNF474", "WDR38", "BASP1", "ZC2HC1A", "TMEM190", "TMEM145", "WDR31")

scc<-all_scores2 %>%
  mutate(phylogenetic_scores = ifelse(phylogenetic_scores > 0.8, "Yes", "No"),
         single_cell_score = ifelse(single_cell_score >= 0.08, "Yes", "No")) %>%
  rename("scRNA-seq" = "single_cell_score", `Comparative Genomics` = "phylogenetic_scores") %>%
  mutate(Type = ifelse(Gene_name %in% ciliaryGenes$Gene_name, "Ciliary", 
                       ifelse(Gene_name %in% negativeCiliary$Gene_name, "Negative", "Candidate"))) %>%
  filter(Gene_name %in% genelist) %>%
  dplyr::select(1,3,6,9) %>%
  rename("Gene name" = "Gene_name")
  
set.seed(30)
scc2<-all_scores2 %>%
  mutate(phylogenetic_scores = ifelse(phylogenetic_scores > 0.8, "Yes", "No"),
         single_cell_score = ifelse(single_cell_score >= 0.08, "Yes", "No")) %>%
  rename("scRNA-seq" = "single_cell_score", `Comparative Genomics` = "phylogenetic_scores") %>%
  mutate(Type = ifelse(Gene_name %in% ciliaryGenes$Gene_name, "Ciliary", 
                       ifelse(Gene_name %in% negativeCiliary$Gene_name, "Negative",
                              ifelse(total_score_sums >= 1, "Candidate", "Unknown")))) %>%
  filter(`Comparative Genomics` == "Yes" & `scRNA-seq` == "No") %>%
  dplyr::select(1,3,6, 9) %>%
  sample_frac(1L) %>%
  slice(1:7) %>%
  rename("Gene name" = "Gene_name")
head(scc2)

scc_all<-rbind(scc, scc2)

scc_all %>% 
  formattable(align =c("l","c","l","c","c","c"),
              list(`Gene Name` = formatter(
                "span", style = ~ style(color = "grey",font.weight = "bold")),
                `Comparative Genomics` = formatter("span",
                                                   style = x ~ style(color = ifelse(x == "Yes", "green", "red")),
                                                   x ~ icontext(ifelse(x == "Yes", "ok", "remove"), ifelse(x == "Yes", "Yes", "No"))),
                `scRNA-seq` = formatter("span",
                                                   style = x ~ style(color = ifelse(x == "Yes", "green", "red")),
                                                   x ~ icontext(ifelse(x == "Yes", "ok", "remove"), ifelse(x == "Yes", "Yes", "No"))),
                `Type` = formatter("span",
                                   style = x ~ style(color = ifelse(x == "Ciliary", "green", ifelse(x == "Negative", "red", "blue"))),
                                   x ~ icontext(ifelse(x == "Ciliary", "ok", ifelse(x == "Negative", "remove", "thumbs-up")), 
                                                ifelse(x == "Ciliary", "Ciliary", ifelse(x == "Negative", "Negative", "Candidate"))))
              ), format = "f", digits = 3) -> tbl

w <- as.htmlwidget(tbl)
w1<-w %>% htmlwidgets::prependContent(htmltools::tags$style("span.highlight {color:red;font-weight: 700;}"))

saveWidget(w1, "table_sc_cg_comparison.html", selfcontained = TRUE)



# 15.11.2022

aday_genes<-all_scores2$Gene_name[all_scores2$total_score_sums > 1]
reyfman_aday<-reyfman[reyfman$Gene_name %in% aday_genes,] # aday = 1082, toplam = 1802, ciliary = 240/688
carraro_aday<-carraro[carraro$Gene_name %in% aday_genes,] # aday = 719, toplam = 1157, ciliary = 231/688
habermann_aday<-habermann[habermann$Gene_name %in% aday_genes,] # aday = 1017, toplam = 2165, ciliary = 287/688
celegans_aday<-celegans[celegans$Gene_name %in% aday_genes,] # aday = 263, toplam = 594, ciliary = 91/688

length(reyfman$Gene_name[reyfman$Gene_name %in% ciliaryGenes$Gene_name])
length(carraro$Gene_name[carraro$Gene_name %in% ciliaryGenes$Gene_name])
length(habermann$Gene_name[habermann$Gene_name %in% ciliaryGenes$Gene_name])
length(celegans$Gene_name[celegans$Gene_name %in% ciliaryGenes$Gene_name])





