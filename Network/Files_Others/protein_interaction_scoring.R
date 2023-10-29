

iref<-fread("All.mitab.06-11-2021.txt", select = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,19,20,29))

xz<-count(iref,taxa)



library(readxl)
library(data.table)
library(dplyr)
library(pbapply)
library(stringr)
library(tidyr)
library(geneName)

source("functions.R")


# Orthology files ----

orthology<-fread("ORTHOLOGY-ALLIANCE_COMBINED_0.tsv.gz") %>% filter(IsBestScore == "Yes")
orthology<-orthology[orthology$Gene2SpeciesName == "Homo sapiens", c(2,4,6,8)]
orthology2<-hgncConverter(orthology, "Gene2Symbol")


#ciliaryGenes<-fread("ciliary_gene_list_updated.txt")
ciliaryGenes<-fread("goldstandart_v2.txt")
ciliaryGenes$scores<-1
ciliaryGenes<-hgncConverter(ciliaryGenes, "Gene_name")
ciliaryGenes$Gene_name[ciliaryGenes$Gene_name == "CILK"]<-"CILK1"
ciliaryGenes$Gene_name[ciliaryGenes$Gene_name == "GPBAR"]<-"GPBAR1"

negativeCiliary<-read_xls("Nevers_2017_NegativeGenesInsightsCiiaryGenes_SuppTable3.xls")
negativeCiliary<-negativeCiliary[,1]
colnames(negativeCiliary)<-"Gene_name"
negativeCiliary<-hgncConverter(negativeCiliary, "Gene_name")

# write.table(ciliaryGenes, "ciliary_gene_list_2022.txt", sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(negativeCiliary, "non_ciliary_gene_list_2022.txt", sep = "\t", row.names = FALSE, quote = FALSE)

ciliatest<-read_xlsx("cilia_genes_references.xlsx")
colnames(ciliatest)[1]<-"Gene_name"
ciliatest<-hgncConverter(ciliatest, "Gene_name")
ciliatest<-anti_join(ciliatest, ciliaryGenes, by = "Gene_name")
ciliatest<-ciliatest %>% dplyr::select(Gene_name) %>% mutate(score = 1)


mouse_orthology<-data.table::fread("HGNC_AllianceHomology.rpt")
colnames(mouse_orthology)[1:16]<-colnames(mouse_orthology)[2:17]
mouse_orthology<-mouse_orthology[,-17]
colnames(mouse_orthology)[2]<-"Gene_name"

mouse_synonyms<-data.table::fread("MRK_List2.rpt")
colnames(mouse_synonyms)[c(7,12)]<-c("Gene_name","Gene_synonyms")
mouse_synonyms<-tidyr::separate_rows(mouse_synonyms, Gene_synonyms, sep = "\\|") %>% 
  dplyr::select(7,12) %>% unique()
mouse_synonyms<-mouse_synonyms[mouse_synonyms$Gene_synonyms != "",]
mouse_orthology<-left_join(mouse_orthology, mouse_synonyms, by = "Gene_name")

hgnc_names<-fread("hgnc_names.txt") %>% filter(Status == "Approved") %>% dplyr::select(1,3) %>% hgncConverter("Approved symbol")
mouse_orthology2<-tidyr::separate_rows(mouse_orthology, `HGNC ID`, sep = "\\|")
mouse_orthology2<-mouse_orthology2[,c(1,2,5,10,16,17)]
mouse_orthology2<-left_join(mouse_orthology2, hgnc_names, by = "HGNC ID")


cele_orthology<-orthology2 %>% filter(Gene1SpeciesName == "Caenorhabditis elegans")


drosophila<-fread("dmel_human_orthologs_disease_fb_2021_06.tsv.gz")
ort_dro<-orthology[orthology$Gene1SpeciesName == "Drosophila melanogaster",]
ort_dro$Gene1Symbol<-gsub("α","alpha",ort_dro$Gene1Symbol)
onlyortdro<-anti_join(ort_dro, drosophila, by = c("Gene1Symbol" = "Dmel_gene_symbol"))

drosophila<-drosophila[,c(2,5)]
colnames(drosophila)<-c("Gene1Symbol","Gene2Symbol")
drosophila<-rbind(drosophila, onlyortdro[,c(1,3)]) %>% unique() %>% hgncConverter("Gene2Symbol")



# Protein - Protein interactions ----


biogrid_mitab<-fread("BIOGRID-ALL-4.4.205.mitab.txt")
biogrid_mitab$`Confidence Values`<-as.numeric(gsub("score:","",biogrid_mitab$`Confidence Values`))

#biogrid_mitab_pro<-biogrid_mitab[grepl("MI:0407|MI:0915|MI:0914|MI:0403", biogrid_mitab$`Interaction Types`),]
biogrid_mitab_pro<-biogrid_mitab

biogrid_mitab_pro_hum<-biogrid_mitab_pro %>%
  filter(`Taxid Interactor A` == "taxid:9606" & `Taxid Interactor B` == "taxid:9606") %>%
  mutate(Interactor_A = str_extract(`Alt IDs Interactor A`, "(?<=locuslink:).*?(?=\\||$)"),
         Interactor_B = str_extract(`Alt IDs Interactor B`, "(?<=locuslink:).*?(?=\\||$)")) %>%
  dplyr::select(c(-1,-2,-3,-4,-5,-6)) %>% 
  hgncConverter("Interactor_B") %>%
  hgncConverter("Interactor_A") %>%
  distinct() %>%
  mutate(Organism = "Homo sapiens")


biogrid_mitab_pro_mus<-biogrid_mitab_pro %>%
  filter(`Taxid Interactor A` == "taxid:10090" & `Taxid Interactor B` == "taxid:10090") %>%
  mutate(Interactor_A = str_extract(`Alt IDs Interactor A`, "(?<=locuslink:).*?(?=\\|entrez|\\||$)"),
         Interactor_B = str_extract(`Alt IDs Interactor B`, "(?<=locuslink:).*?(?=\\|entrez|\\||$)")) %>%
  dplyr::select(c(-1,-2,-3,-4,-5,-6)) %>%
  distinct()

biogrid_mitab_pro_mus$Interactor_A1<-mouse_orthology2$`Approved symbol`[match(biogrid_mitab_pro_mus$Interactor_A,
                                                                             mouse_orthology2$Gene_name)]
biogrid_mitab_pro_mus$Interactor_A2<-mouse_orthology2$`Approved symbol`[match(biogrid_mitab_pro_mus$Interactor_A,
                                                                              mouse_orthology2$Gene_synonyms)]
biogrid_mitab_pro_mus$Interactor_A<-ifelse(!is.na(biogrid_mitab_pro_mus$Interactor_A1), 
                                           biogrid_mitab_pro_mus$Interactor_A1, 
                                           biogrid_mitab_pro_mus$Interactor_A2)
biogrid_mitab_pro_mus<-biogrid_mitab_pro_mus[,c(-12,-13)]

biogrid_mitab_pro_mus$Interactor_B1<-mouse_orthology2$`Approved symbol`[match(biogrid_mitab_pro_mus$Interactor_B,
                                                                              mouse_orthology2$Gene_name)]
biogrid_mitab_pro_mus$Interactor_B2<-mouse_orthology2$`Approved symbol`[match(biogrid_mitab_pro_mus$Interactor_B,
                                                                              mouse_orthology2$Gene_synonyms)]
biogrid_mitab_pro_mus$Interactor_B<-ifelse(!is.na(biogrid_mitab_pro_mus$Interactor_B1), 
                                           biogrid_mitab_pro_mus$Interactor_B1, 
                                           biogrid_mitab_pro_mus$Interactor_B2)
biogrid_mitab_pro_mus<-biogrid_mitab_pro_mus[,c(-12,-13)]
biogrid_mitab_pro_mus$Organism<-"Mus musculus"



biogrid_mitab_pro_dro<-biogrid_mitab_pro %>%
  filter(`Taxid Interactor A` == "taxid:7227" & `Taxid Interactor B` == "taxid:7227") %>%
  mutate(Interactor_A = str_extract(`Alt IDs Interactor A`, "(?<=locuslink:).*?(?=\\|entrez|\\||$)"),
         Interactor_B = str_extract(`Alt IDs Interactor B`, "(?<=locuslink:).*?(?=\\|entrez|\\||$)")) %>%
  dplyr::select(c(-1,-2,-3,-4,-5,-6)) %>%
  distinct()

biogrid_mitab_pro_dro$Interactor_A<-gsub('"','',biogrid_mitab_pro_dro$Interactor_A)
biogrid_mitab_pro_dro$Interactor_B<-gsub('"','',biogrid_mitab_pro_dro$Interactor_B)

biogrid_mitab_pro_dro$Interactor_A<-drosophila$Gene2Symbol[match(biogrid_mitab_pro_dro$Interactor_A,
                                                                              drosophila$Gene1Symbol)]

biogrid_mitab_pro_dro$Interactor_B<-drosophila$Gene2Symbol[match(biogrid_mitab_pro_dro$Interactor_B,
                                                                  drosophila$Gene1Symbol)]
biogrid_mitab_pro_dro$Organism<-"Drosophila melanogaster"


biogrid_mitab_pro_cele<-biogrid_mitab_pro %>%
  filter(`Taxid Interactor A` == "taxid:6239" & `Taxid Interactor B` == "taxid:6239") %>%
  mutate(Interactor_A = str_extract(`Alt IDs Interactor A`, "(?<=locuslink:).*?(?=\\|entrez|\\||$)"),
         Interactor_B = str_extract(`Alt IDs Interactor B`, "(?<=locuslink:).*?(?=\\|entrez|\\||$)")) %>%
  dplyr::select(c(-1,-2,-3,-4,-5,-6)) %>%
  distinct()

biogrid_mitab_pro_cele$Interactor_A<-gsub("CELE_","",biogrid_mitab_pro_cele$Interactor_A)
biogrid_mitab_pro_cele$Interactor_B<-gsub("CELE_","",biogrid_mitab_pro_cele$Interactor_B)


biogrid_mitab_pro_cele$Interactor_A<-cele_orthology$Gene2Symbol[match(biogrid_mitab_pro_cele$Interactor_A,
                                                                      cele_orthology$Gene1Symbol)]

biogrid_mitab_pro_cele$Interactor_B<-cele_orthology$Gene2Symbol[match(biogrid_mitab_pro_cele$Interactor_B,
                                                                      cele_orthology$Gene1Symbol)]
biogrid_mitab_pro_cele$Organism<-"Caenorhabditis elegans"


mit_biogrid<-rbind(biogrid_mitab_pro_hum, biogrid_mitab_pro_mus, biogrid_mitab_pro_dro, biogrid_mitab_pro_cele)
mit_biogrid<-mit_biogrid[!is.na(mit_biogrid$Interactor_A) & !is.na(mit_biogrid$Interactor_B)]



 


# Intact ----

intact<-fread("intact.txt.1", select = c(5,6,7,8,9,10,11,12,13,15,21,22)) %>% unique() %>%
  mutate(ppi_score = as.numeric(str_extract(`Confidence value(s)`, "(?<=intact-miscore:).+"))) %>%
  filter(!is.na(ppi_score)) %>%
  mutate(`Publication Identifier(s)` = str_extract(`Publication Identifier(s)`, "pubmed:.*?(?=\\||$)"))

#intact_pro<-intact[grepl("MI:0407|MI:0915|MI:0914|MI:0403", intact$`Interaction type(s)`),]
intact_pro<-intact
intact_hum<-intact_pro %>% filter(`Taxid interactor A` == "taxid:9606(human)|taxid:9606(Homo sapiens)",
                              `Taxid interactor B` == "taxid:9606(human)|taxid:9606(Homo sapiens)",
                              grepl("gene name", `Alias(es) interactor A`),
                              grepl("gene name", `Alias(es) interactor B`))
intact_hum1<-intact_hum %>%
  mutate(`Alias(es) interactor A` = gsub(".*uniprotkb:(.+)\\(gene name).*", "\\1", intact_hum$`Alias(es) interactor A`)) %>%
  mutate(`Alias(es) interactor B` = gsub(".*uniprotkb:(.+)\\(gene name).*", "\\1", intact_hum$`Alias(es) interactor B`))

intact_hum1<-hgncConverter(intact_hum1, "Alias(es) interactor B")
intact_hum1<-hgncConverter(intact_hum1, "Alias(es) interactor A")

intact_hum1$Organism<-"Homo sapiens"

## Intact celegans ----

intact_cele<-intact_pro %>% filter(`Taxid interactor A` == "taxid:6239(caeel)|taxid:6239(Caenorhabditis elegans)",
                               `Taxid interactor B` == "taxid:6239(caeel)|taxid:6239(Caenorhabditis elegans)",
                               grepl("orf name", `Alias(es) interactor A`) | grepl("gene name", `Alias(es) interactor A`),
                               grepl("orf name", `Alias(es) interactor B`) | grepl("gene name", `Alias(es) interactor B`))

intact_cele[["Alias(es) interactor A"]]<-ifelse(grepl("\\(gene name\\)", intact_cele$`Alias(es) interactor A`),
                                                gsub(".*uniprotkb:(.+)\\(gene name).*", "\\1", intact_cele$`Alias(es) interactor A`),
                                                ifelse(grepl("CELE", intact_cele$`Alias(es) interactor A`),
                                                       str_extract(intact_cele$`Alias(es) interactor A`, "(?<=CELE_).+?(?=\\(orf name\\))"),
                                                       str_extract(intact_cele$`Alias(es) interactor A`, "(?<=uniprotkb:)(?!.*synonym).+?(?=\\(orf name\\))")))


intact_cele[["Alias(es) interactor B"]]<-ifelse(grepl("\\(gene name\\)", intact_cele$`Alias(es) interactor B`),
                                                gsub(".*uniprotkb:(.+)\\(gene name).*", "\\1", intact_cele$`Alias(es) interactor B`),
                                                ifelse(grepl("CELE", intact_cele$`Alias(es) interactor B`),
                                                       str_extract(intact_cele$`Alias(es) interactor B`, "(?<=CELE_).+?(?=\\(orf name\\))"),
                                                       str_extract(intact_cele$`Alias(es) interactor B`, "(?<=uniprotkb:)(?!.*synonym).+?(?=\\(orf name\\))")))

intact_cele1<-intact_cele
intact_cele1$`Alias(es) interactor A`<-cele_orthology$Gene2Symbol[match(intact_cele1$`Alias(es) interactor A`, cele_orthology$Gene1Symbol)]
intact_cele1$`Alias(es) interactor B`<-cele_orthology$Gene2Symbol[match(intact_cele1$`Alias(es) interactor B`, cele_orthology$Gene1Symbol)]
intact_cele1<-na.omit(intact_cele1)
intact_cele1$Organism<-"Caenorhabditis elegans"



## Intact mouse ----

intact_mus<-intact %>% filter(`Taxid interactor A` == "taxid:10090(mouse)|taxid:10090(Mus musculus)",
                              `Taxid interactor B` == "taxid:10090(mouse)|taxid:10090(Mus musculus)")

intact_mus1<-intact_mus %>% 
  mutate(`Alias(es) interactor A` = gsub(".*uniprotkb:(.+)\\(gene name).*", "\\1", intact_mus$`Alias(es) interactor A`)) %>%
  mutate(`Alias(es) interactor B` = gsub(".*uniprotkb:(.+)\\(gene name).*", "\\1", intact_mus$`Alias(es) interactor B`))

intact_mus1<-intact_mus1 %>%
  mutate(`Alias(es) interactor A` = ifelse(nchar(intact_mus1$`Alias(es) interactor A`)<15, `Alias(es) interactor A`,
                                           str_extract(`Alias(es) interactor A`, "(?<=psi-mi:).+?(?=_mouse_gene)"))) %>%
  mutate(`Alias(es) interactor B` = ifelse(nchar(intact_mus1$`Alias(es) interactor B`)<15, `Alias(es) interactor B`,
                                           str_extract(`Alias(es) interactor B`, "(?<=psi-mi:).+?(?=_mouse_gene)")))

intact_mus1<-intact_mus1[!is.na(intact_mus1$`Alias(es) interactor A`) | !is.na(intact_mus1$`Alias(es) interactor A`),]
intact_mus1$`Alias(es) interactor A`<-str_to_title(intact_mus1$`Alias(es) interactor A`)
intact_mus1$`Alias(es) interactor B`<-str_to_title(intact_mus1$`Alias(es) interactor B`)

intact_mus1$`Alias(es) interactor A`<-mouse_orthology2$`Approved symbol`[match(intact_mus1$`Alias(es) interactor A`, mouse_orthology2$Gene_name)]
intact_mus1$`Alias(es) interactor B`<-mouse_orthology2$`Approved symbol`[match(intact_mus1$`Alias(es) interactor B`, mouse_orthology2$Gene_name)]
intact_mus1<-na.omit(intact_mus1)
intact_mus1$Organism<-"Mus musculus"


## Intact drosophila ----

intact_dro<-intact %>% 
  filter(`Taxid interactor A` == 'taxid:7227(drome)|taxid:7227("Drosophila melanogaster (Fruit fly)")') %>%
  filter(`Taxid interactor B` == 'taxid:7227(drome)|taxid:7227("Drosophila melanogaster (Fruit fly)")')

intact_dro[["Alias(es) interactor A"]]<-gsub(".*uniprotkb:(.+)\\(gene name).*", "\\1", intact_dro$`Alias(es) interactor A`)
intact_dro[["Alias(es) interactor A"]]<-gsub("Dmel\\\\","",intact_dro$`Alias(es) interactor A`)
intact_dro[["Alias(es) interactor A"]]<-gsub("BEST:","",intact_dro$`Alias(es) interactor A`)
intact_dro[["Alias(es) interactor A"]]<-gsub('"',"",intact_dro$`Alias(es) interactor A`)

intact_dro[["Alias(es) interactor B"]]<-gsub(".*uniprotkb:(.+)\\(gene name).*", "\\1", intact_dro$`Alias(es) interactor B`)
intact_dro[["Alias(es) interactor B"]]<-gsub("Dmel\\\\","",intact_dro$`Alias(es) interactor B`)
intact_dro[["Alias(es) interactor B"]]<-gsub("BEST:","",intact_dro$`Alias(es) interactor B`)
intact_dro[["Alias(es) interactor B"]]<-gsub('"',"",intact_dro$`Alias(es) interactor B`)


intact_dro1<-intact_dro
intact_dro1$`Alias(es) interactor A`<-drosophila$Gene2Symbol[match(intact_dro1$`Alias(es) interactor A`, drosophila$Gene1Symbol)]
intact_dro1$`Alias(es) interactor B`<-drosophila$Gene2Symbol[match(intact_dro1$`Alias(es) interactor B`, drosophila$Gene1Symbol)]
intact_dro1<-na.omit(intact_dro1)
intact_dro1$Organism<-"Drosophila melanogaster"



intact_all<-rbind(intact_hum1, intact_mus1, intact_dro1, intact_cele1)

onlyintact<-anti_join(intact_all, mit_biogrid, by = c("Alias(es) interactor A" = "Interactor_A", 
                                                "Alias(es) interactor B" = "Interactor_B",
                                                "Publication Identifier(s)" = "Publication Identifiers"))

onlyintact<-onlyintact[,c(-10,-11,-12)]
mit_biogrid<-mit_biogrid[,-10]
colnames(onlyintact)<-colnames(mit_biogrid)

bioin<-rbind(mit_biogrid, onlyintact)


## Huri ----

library(biomaRt)

huri<-fread("HI-union.psi", select = c(3,4,7,8,9,10,11,12,15), header = FALSE) %>%
  mutate(V3 = str_extract(V3, "ENSG.*?(?=\\.)"),
         V4 = str_extract(V4, "ENSG.*?(?=\\.)"),
         V15 = as.numeric(str_extract(V15, "(?<=author score:).*")),
         `Source Database` = "HuRI",
         Organism = "Homo sapiens")

huri$V9[grepl("unassigned", huri$V9)]<-"pubmed:32296183"


ensembl<-useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
mart<-useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

hi_genes<-getBM(attributes = c("ensembl_gene_id","external_gene_name"), 
                filters = "ensembl_gene_id",
                values = unique(c(huri$V3, huri$V4)),
                mart = mart) %>% na.omit()

hi_genes<-hgncConverter(hi_genes, "external_gene_name")
hi_genes$external_gene_name[hi_genes$external_gene_name == ""]<-hi_genes$ensembl_gene_id[hi_genes$external_gene_name == ""]

huri$V3<-hi_genes$external_gene_name[match(huri$V3, hi_genes$ensembl_gene_id)]
huri$V4<-hi_genes$external_gene_name[match(huri$V4, hi_genes$ensembl_gene_id)]

colnames(huri)[1:8]<-colnames(bioin)[1:8]
colnames(huri)[9]<-"Confidence Values"

onlyhuri<-anti_join(huri, bioin, by = c("Interactor_A", "Interactor_B", "Publication Identifiers"))

bioinhu<-rbind(bioin, onlyhuri)
#bioinhu_pro<-bioinhu
bioinhu_pro<-bioinhu[!grepl("MI:0254", bioinhu$`Interaction Detection Method`),]




library(ontologyIndex)

owl<-get_ontology("mi.owl")


idmpar<-data.frame(method = c("MI:0013","MI:0090","MI:0254","MI:0255","MI:0401","MI:0428"),
                   weight = c(1,0.66,0.10,0.10,1,0.33))

itpar<-data.frame(method = c("MI:0208","MI:0403","MI:0914","MI:0915","MI:0407"),
                  weight = c(0.10,0.33,0.33,0.66,1),
                  group = c(3,1,2,2,2))

# idm<-data.frame(method = unique(biogrid_mitab_pro$`Interaction Detection Method`),
#                 weight = c(0.66,1,1,1,1,1,1,0.05,0.33,0.66,1,1))
# 
# it<-data.frame(method = unique(biogrid_mitab_pro$`Interaction Types`),
#                weight = c(1,0.66,0.33,0.33))


# Genetic interactions ----

bioinhu_gen<-bioinhu[grepl("MI:0254", bioinhu$`Interaction Detection Method`),-10]


alliance_dro<-fread("dro_alliance.tsv.gz", select = c(5:13))
alliance_dro$Interactor_A<-str_extract(alliance_dro$`Alias(es) interactor A`, "(?<=flybase:).+?(?=\\(gene name\\))")
alliance_dro$Interactor_A<-gsub('"',"",alliance_dro$Interactor_A)
alliance_dro$Interactor_B<-str_extract(alliance_dro$`Alias(es) interactor B`, "(?<=flybase:).+?(?=\\(gene name\\))")
alliance_dro$Interactor_B<-gsub('"',"",alliance_dro$Interactor_B)
alliance_dro$Interactor_A<-drosophila$Gene2Symbol[match(alliance_dro$Interactor_A, drosophila$Gene1Symbol)]
alliance_dro$Interactor_B<-drosophila$Gene2Symbol[match(alliance_dro$Interactor_B, drosophila$Gene1Symbol)]
alliance_dro<-alliance_dro[!(is.na(alliance_dro$Interactor_A) | is.na(alliance_dro$Interactor_B)),]

alliance_dro$`Publication Identifier(s)`<-str_extract(alliance_dro$`Publication Identifier(s)`, "pubmed:.+?(?=$)")
alliance_dro$`Publication Identifier(s)`[is.na(alliance_dro$`Publication Identifier(s)`)]<-"Unknown"


onlyalldro<-anti_join(alliance_dro, bioinhu_gen, 
                      by = c("Publication Identifier(s)" = "Publication Identifiers", 
                             "Interactor_A", "Interactor_B")) %>%
  dplyr::select(c(10,11,3:9)) %>%
  mutate(Organism = "Drosophila melanogaster")

colnames(onlyalldro)<-colnames(bioinhu_gen)  


alliance_cele<-fread("cele_alliance.tsv.gz", select = c(5:13))

alliance_cele$Interactor_A<-str_extract(alliance_cele$`Alias(es) interactor A`, "(?<=wormbase:).+?(?=\\(public_name\\))")
alliance_cele$Interactor_B<-str_extract(alliance_cele$`Alias(es) interactor B`, "(?<=wormbase:).+?(?=\\(public_name\\))")
alliance_cele$Interactor_A<-cele_orthology$Gene2Symbol[match(alliance_cele$Interactor_A, cele_orthology$Gene1Symbol)]
alliance_cele$Interactor_B<-cele_orthology$Gene2Symbol[match(alliance_cele$Interactor_B, cele_orthology$Gene1Symbol)]
alliance_cele<-alliance_cele[!(is.na(alliance_cele$Interactor_A) | is.na(alliance_cele$Interactor_B)),]

onlyallcele<-anti_join(alliance_cele, bioinhu_gen, 
                       by = c("Publication Identifier(s)" = "Publication Identifiers", 
                              "Interactor_A", "Interactor_B")) %>%
  dplyr::select(c(10,11,3:9)) %>%
  mutate(Organism = "Caenorhabditis elegans")

colnames(onlyallcele)<-colnames(bioinhu_gen)  


genetic_int_all<-rbind(bioinhu_gen, onlyalldro, onlyallcele)


bioinhu_pro<-rbind(bioinhu_pro[,-10], genetic_int_all)

bioinhu1<-bioinhu_pro %>% 
  mutate(method = str_extract(`Interaction Detection Method`, '(?<=psi-mi:").*?(?=")'),
         type = str_extract(`Interaction Types`, '(?<=psi-mi:").*?(?=")'))



allmethods<-owl$ancestors[bioinhu1$method]
methodlist<-lapply(allmethods, function(x) idmpar$weight[which(idmpar$method %in% x)])
is.na(methodlist)<-lengths(methodlist) == 0
methodlist<-sapply(methodlist, "[[", 1)
methodlist<-unlist(methodlist, use.names = FALSE)


alltypes<-owl$ancestors[bioinhu1$type]
typelist<-lapply(alltypes, function(x) itpar$weight[which(itpar$method %in% x)])
is.na(typelist)<-lengths(typelist) == 0
typelist<-lapply(typelist, max)
typelist<-sapply(typelist, "[[", 1)
typelist<-unlist(typelist, use.names = FALSE)

grouplist<-lapply(alltypes, function(x) itpar$group[which(itpar$method %in% x)])
is.na(grouplist)<-lengths(grouplist) == 0
grouplist<-lapply(grouplist, max)
grouplist<-sapply(grouplist, "[[", 1)
grouplist<-unlist(grouplist, use.names = FALSE)


bioinhu1$idm<-methodlist
bioinhu1$it<-typelist
bioinhu1$it2<-grouplist
#bioinhu1<-bioinhu1[!is.na(bioinhu1$it),]
#bioinhu1<-bioinhu1[!is.na(bioinhu1$idm),]
bioinhu1$it[is.na(bioinhu1$it)]<-0.05
bioinhu1$idm[is.na(bioinhu1$idm)]<-0.05
bioinhu1$it2[is.na(bioinhu1$it2)]<-3
bioinhu1<-unique(bioinhu1)

# bioinhu2<-bioinhu1 %>% distinct(Interactor_A, Interactor_B, `Publication Identifiers`, 
#                                 `Interaction Detection Method`, `Interaction Types`, .keep_all = TRUE)
# 
# bioinhu3<-bioinhu1 %>% distinct(Interactor_A, Interactor_B, `Publication Identifiers`, 
#                                 `Interaction Detection Method`, .keep_all = TRUE)
# 
# bioinhu4<-bioinhu1 %>% distinct(Interactor_A, Interactor_B, `Publication Identifiers`, 
#                                 `Interaction Types`, .keep_all = TRUE)

# yy<-bioinhu1 %>% group_by(Interactor_A, Interactor_B, `Publication Identifiers`) %>%
#   add_count(`Interaction Detection Method`) %>%
#   add_count(`Interaction Types`) %>%
#   ungroup() %>%
#   mutate(a = idm * n) %>%
#   mutate(b = it * nn) %>%
#   group_by(Interactor_A, Interactor_B) %>%
#   mutate(a1 = sum(a), gscv = sum(max(idm))) %>%
#   mutate(b1 = sum(b), gscv2 = sum(max(it))) %>%
#   ungroup() %>%
#   mutate(method_score = log(a1 + 1, a1 + gscv + 1)) %>%
#   mutate(type_score = log(b1 + 1, b1 + gscv2 + 1))



yy1<-bioinhu1 %>% 
  group_by(Interactor_A, Interactor_B, `Interaction Detection Method`) %>%
  mutate(n = n_distinct(`Publication Identifiers`)) %>%
  ungroup() %>%
  group_by(Interactor_A, Interactor_B, `Interaction Types`) %>%
  mutate(nn = n_distinct(`Publication Identifiers`)) %>%
  ungroup() %>%
  mutate(a = idm * n, 
         b = it * nn)

yy1<-yy1[!is.na(yy1$Interactor_A),]
yy1<-yy1[!is.na(yy1$Interactor_B),]

yy2<-yy1 %>%
  group_by(Interactor_A, Interactor_B) %>%
  distinct(`Interaction Detection Method`, .keep_all = TRUE) %>%
  mutate(a1 = sum(a))

yy3<-yy1 %>%
  group_by(Interactor_A, Interactor_B) %>%
  distinct(`Interaction Types`, .keep_all = TRUE) %>%
  mutate(b1 = sum(b))

yy4<-yy1 %>% left_join(yy2, by = names(yy1)) %>% left_join(yy3, by = names(yy1)) %>%
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


yy5<-yy4[,c(-11:-23)]

yy5$ppi_score<-(sqrt(yy5$method_score * yy5$type_score))

yy5<-as.data.frame(yy5)
# yy5 ----

ciliaryGenes<-na.omit(ciliaryGenes)

yy6<-yy5 %>%
  left_join(ciliaryGenes, by = c("Interactor_B" = "Gene_name")) %>%
  mutate(across(c("scores"), ~replace_na(.x, 0))) %>%
  rename("Ciliary_A" = "scores")

yy6<-yy6 %>%
  left_join(ciliaryGenes, by = c("Interactor_A" = "Gene_name")) %>%
  mutate(across(c("scores"), ~replace_na(.x, 0))) %>%
  rename("Ciliary_B" = "scores") %>%
  mutate(Ciliary_A = Ciliary_A * ppi_score,
         Ciliary_B = Ciliary_B * ppi_score)



int_a<-unique(yy6[,c(1,2,14)])
int_b<-unique(yy6[,c(2,1,15)])
colnames(int_b)<-colnames(int_a)

int_all<-rbind(int_a, int_b)

alist<-int_all[,c(1,3)]
alist1x<-alist %>% group_by(Interactor_A) %>% 
  mutate(gscv = sum(max(Ciliary_A)), a1 = sum(Ciliary_A), ln = length(Interactor_A)) %>%
  mutate(ln1 = log((ln / 5)+1)) %>%
  #mutate(ln = ifelse(ln<=30, ln, 30)) %>%
  # mutate(score = log((ln1 %/% 5) +1) * (exp(mean(log(Ciliary_A+1)))-1)) %>%
  # mutate(score2 = ln1 %/% 5 * (exp(mean(log(Ciliary_A+1)))-1)) %>%
  mutate(scoreFinal = ln1 * (exp(mean(log(Ciliary_A+1)))-1))
  # mutate(score3 = mean(Ciliary_A)) %>%
  # mutate(score4 = exp(mean(log(Ciliary_A+1)))-1)


library(quantable)


alist2x<-unique(alist1x[,c(1,7)])

colnames(alist2x)<-c("Gene_name","interaction_scores")

alist2x$interaction_scores<-normalization0(log(alist2x$interaction_scores+1)-1)

hgnc_names<-fread("hgnc_names.txt", select = c(1,3))
hgnc_names<-unique(hgnc_names)

prot_int<-left_join(hgnc_names[,2], alist2x, by = c("Approved symbol" = "Gene_name"))
colnames(prot_int)<-c("Gene_name", "interaction_scores")
prot_int$interaction_scores[is.na(prot_int$interaction_scores)]<-0

write.table(prot_int, "protein_interaction_scores_2022.txt", row.names = FALSE, quote = FALSE, sep = "\t")
#abc<-anti_join(alist2x, all_genes, by = "Gene_name")

# ****************************************************************************

#x<-robustscale(alist2x[alist2x$scoreFinal != 0, 5])
#alist2x$norm[alist2x$scoreFinal != 0]<-x$data$scoreFinal

# alist2x$log<-log(alist2x$norm + 1)
# alist2x$log_ln<-log((alist2x$ln %/% 5)+1)
# plot(density(alist2x$log))

cilprot<-semi_join(alist2x, ciliatest, by = c("Interactor_A" = "Gene_name"))
nocilprot<-semi_join(alist2x, negativeCiliary, by = c("Interactor_A" = "Gene_name"))

cilprot$is.ciliary<-1
nocilprot$is.ciliary<-0
cilnoprot<-rbind(cilprot, nocilprot)
cilnoprot$predicted<-cilnoprot$scoreFinal
# cilnoprot$predicted<-ifelse(cilnoprot$scoreFinal >= 0.1, 1, 0)
# 
# cilnoprot$is.ciliary<-as.factor(cilnoprot$is.ciliary)
# cilnoprot$predicted<-as.factor(cilnoprot$predicted)
# 
# cnfs<-confusionMatrix(cilnoprot$predicted, cilnoprot$is.ciliary, positive = "1", mode="prec_recall")
# cnfs$byClass[["F1"]]

classifier(cilnoprot, "scoreFinal", "is.ciliary")



### protein scores ----

alisttest<-semi_join(alist2x, ciliatest, by = c("Interactor_A" = "Gene_name"))
alisttest2<-semi_join(alist2x, ciliaryGenes, by = c("Interactor_A" = "Gene_name"))


alist1<-alist[alist$Ciliary_B>0,] %>% group_by(Interactor_A) %>% mutate(score_a = exp(mean(log(Ciliary_B))))
alist1<-unique(alist1[,c(1,3)])






## Scoring genetic interactions ----


bioflyworm<-genetic_int_all %>% 
  mutate(method = str_extract(`Interaction Detection Method`, '(?<=psi-mi:").*?(?=")'),
         type = str_extract(`Interaction Types`, '(?<=psi-mi:").*?(?=")'))



allmethods<-owl$ancestors[bioflyworm$method]
methodlist<-lapply(allmethods, function(x) idmpar$weight[which(idmpar$method %in% x)])
is.na(methodlist)<-lengths(methodlist) == 0
methodlist<-sapply(methodlist, "[[", 1)
methodlist<-unlist(methodlist, use.names = FALSE)


alltypes<-owl$ancestors[bioflyworm$type]
typelist<-lapply(alltypes, function(x) itpar$weight[which(itpar$method %in% x)])
is.na(typelist)<-lengths(typelist) == 0
typelist<-lapply(typelist, max)
typelist<-sapply(typelist, "[[", 1)
typelist<-unlist(typelist, use.names = FALSE)

grouplist<-lapply(alltypes, function(x) itpar$group[which(itpar$method %in% x)])
is.na(grouplist)<-lengths(grouplist) == 0
grouplist<-lapply(grouplist, max)
grouplist<-sapply(grouplist, "[[", 1)
grouplist<-unlist(grouplist, use.names = FALSE)


bioflyworm$idm<-methodlist
bioflyworm$it<-typelist
bioflyworm$it2<-grouplist
#bioflyworm<-bioflyworm[!is.na(bioflyworm$it),]
#bioflyworm<-bioflyworm[!is.na(bioflyworm$idm),]
bioflyworm$it[is.na(bioflyworm$it)]<-0.05
bioflyworm$idm[is.na(bioflyworm$idm)]<-0.05
bioflyworm$it2[is.na(bioflyworm$it2)]<-4
bioflyworm<-unique(bioflyworm)


xx1<-bioflyworm %>% 
  group_by(Interactor_A, Interactor_B, `Interaction Detection Method`) %>%
  mutate(n = n_distinct(`Publication Identifiers`)) %>%
  ungroup() %>%
  group_by(Interactor_A, Interactor_B, `Interaction Types`) %>%
  mutate(nn = n_distinct(`Publication Identifiers`)) %>%
  ungroup() %>%
  mutate(a = idm * n, 
         b = it * nn)

xx1<-xx1[!is.na(xx1$Interactor_A),]

xx2<-xx1 %>%
  group_by(Interactor_A, Interactor_B) %>%
  distinct(`Interaction Detection Method`, .keep_all = TRUE) %>%
  mutate(a1 = sum(a))

xx3<-xx1 %>%
  group_by(Interactor_A, Interactor_B) %>%
  distinct(`Interaction Types`, .keep_all = TRUE) %>%
  mutate(b1 = sum(b))

xx4<-xx1 %>% left_join(xx2, by = names(xx1)) %>% left_join(xx3, by = names(xx1)) %>%
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


xx5<-xx4[,c(-11:-23)]

xx5$ppi_score<-(sqrt(xx5$method_score * xx5$type_score))


xx6<-xx5 %>%
  left_join(ciliaryGenes, by = c("Interactor_A" = "Gene_name")) %>%
  left_join(ciliaryGenes, by = c("Interactor_B" = "Gene_name")) %>%
  rename("Ciliary_A" = "scores.x", "Ciliary_B" = "scores.y") %>%
  mutate(across(c("Ciliary_A", "Ciliary_B"), ~replace_na(.x, 0))) %>%
  mutate(Ciliary_A = Ciliary_A * ppi_score,
         Ciliary_B = Ciliary_B * ppi_score)


int_a<-unique(xx6[,c(1,2,15)])
int_b<-unique(xx6[,c(2,1,14)])
colnames(int_b)<-colnames(int_a)

int_all<-rbind(int_a, int_b)

blist<-int_all[,c(1,3)]
blist1x<-blist %>% group_by(Interactor_A) %>% 
  mutate(gscv = sum(max(Ciliary_B)), a1 = sum(Ciliary_B), ln = length(Interactor_A)/5) %>%
  mutate(ln = ifelse(ln<=30, ln, 30)) %>%
  mutate(score = ln %/% 5 * (exp(mean(log(Ciliary_B+1)))-1)) %>%
  ungroup() #%>%


blist2x<-unique(blist1x[,c(1,6)])


blisttest<-semi_join(blist2x, ciliatest, by = c("Interactor_A" = "Gene_name"))
blisttest2<-semi_join(blist2x, ciliaryGenes, by = c("Interactor_A" = "Gene_name"))


blist1<-alist[alist$Ciliary_B>0,] %>% group_by(Interactor_A) %>% mutate(score_a = exp(mean(log(Ciliary_B))))
blist1<-unique(alist1[,c(1,3)])


