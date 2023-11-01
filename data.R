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


