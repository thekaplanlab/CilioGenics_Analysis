

library(data.table)
library(XML)
library(dplyr)
library(clipr)
library(xml2)

hgnc2<-fread("hgnc_custom.txt") %>%
  filter(Status == "Approved") %>%
  dplyr::select(c(-1,-2))
a<-hgnc2[,1:2]
b<-hgnc2[,c(1,3)]
colnames(b)<-colnames(a)
hgnc2<-unique(rbind(a,b))
colnames(hgnc2)<-c("Gene_name", "Gene_synonyms")


protlas<-xmlParse("proteinatlas.xml.gz")

data<-read_xml("http://www.proteinatlas.org/ENSG00000000971.xml")
data_xml<-xmlParse(data)

somefun = function(x, ...) {
  #This is a "branch" function
  #x is a XML node - everything inside element <MedlineCitation>
  # find element <ArticleTitle> inside and print it:
  ns <- getNodeSet(x,path = "//entry")
  value <- xmlValue(ns[[1]])
  print(value)
}

aa<-xmlEventParse(
  file = "http://www.proteinatlas.org/ENSG00000097007.xml", 
  handlers = NULL, 
  branches = list(entry = somefun)
)
edge.nodes <- xml_find_all( data, ".//cellType")
dat<-xmlToDataFrame(nodes=getNodeSet(data_xml, "//tissueCell"))

dat<-getNodeSet(data_xml, "//cellType")
xmlParent(data_xml, "//cellType")

proteinatlas<-fread("proteinatlas.tsv") %>%
  dplyr::select(c(1,3,4:12))
proatlas1<-as.data.table(proteinatlas)

pb<-pbapply::timerProgressBar(min = 1, max = length(proteinatlas$Ensembl), style = 2)
for (i in 1:length(proteinatlas$Ensembl)){
  ensgene<-proteinatlas$Ensembl[i]
  link<-paste0("http://www.proteinatlas.org/", ensgene, ".xml")
  data <- xmlParse(link)
  xml_data<-xmlToList(data)
  cdata<-xpathSApply(data, "//text()", xmlValue)
  cdata<-gsub("\n","",cdata)
  a<-cdata[grepl("cilia", cdata, ignore.case = TRUE)]
  if (length(a)>0){
    proatlas1[i,1:=proteinatlas$Gene[i]]
    proatlas1[i,2:=paste(a, collapse = ",")]
    proatlas1[i,7:="YES"]
  }
  else {proatlas1[i,2:=NA]
    proatlas1[i,7:="NO"]}
  
  a<-cdata[grepl("cilium", cdata, ignore.case = TRUE)]
  if (length(a)>0){
    proatlas1[i,1:=proteinatlas$Gene[i]]
    proatlas1[i,3:=paste(a, collapse = ",")]
    proatlas1[i,8:="YES"]
  }
  else {proatlas1[i,3:=NA]
    proatlas1[i,8:="NO"]}
  
  a<-cdata[grepl("centrosome", cdata, ignore.case = TRUE)]
  if (length(a)>0){
    proatlas1[i,1:=proteinatlas$Gene[i]]
    proatlas1[i,4:=paste(a, collapse = ",")]
    proatlas1[i,9:="YES"]
  }
  else {proatlas1[i,4:=NA]
    proatlas1[i,9:="NO"]}
  
  a<-cdata[grepl("flagella", cdata, ignore.case = TRUE)]
  if (length(a)>0){
    proatlas1[i,1:=proteinatlas$Gene[i]]
    proatlas1[i,5:=paste(a, collapse = ",")]
    proatlas1[i,10:="YES"]
  }
  else {proatlas1[i,5:=NA]
    proatlas1[i,10:="NO"]}
  
  a<-cdata[grepl("flagellum", cdata, ignore.case = TRUE)]
  if (length(a)>0){
    proatlas1[i,1:=proteinatlas$Gene[i]]
    proatlas1[i,6:=paste(a, collapse = ",")]
    proatlas1[i,11:="YES"]
  }
  else {proatlas1[i,6:=NA]
    proatlas1[i,11:="NO"]}
  
  pbapply::setTimerProgressBar(pb, i)
}
write.table(proatlas1, "protein_atlas_2022.txt", quote = FALSE, sep = "\t", row.names = FALSE)


proatlas2<-fread("protein_atlas_2022.txt")
colnames(proatlas2)<-c("Gene_name", "Cilia_comments", "Cilium_comments", "Centrosome_comments", "Flagella_comments", "Flagellum_comments",
                       "Cilia", "Cilium", "Centrosome", "Flagella", "Flagellum")

proatlas3<- subset(proatlas2, Cilia == "YES" | Cilium == "YES" | Centrosome == "YES" | Flagella == "YES" | Flagellum == "YES")

write_clip(proatlas3)


progenes<-proatlas3[,1]

progenes$Centrosome<-ifelse(apply(proatlas3[,2:6], 1, function(x) any(grepl("centrosome", x, ignore.case = TRUE))), "YES", "NO")
progenes$Cilia_expression<-ifelse(apply(proatlas3[,2:6], 1, function(x) any(grepl("cilia|cilium", x, ignore.case = TRUE))), "YES", "NO")
progenes$Flagella<-ifelse(apply(proatlas3[,2:6], 1, function(x) any(grepl("flagella|flagellum", x, ignore.case = TRUE))), "YES", "NO")
progenes$Cilia_localization<-ifelse(apply(proatlas3[,2:6], 1, function(x) any(grepl("cilia|cilium", x, ignore.case = TRUE) & grepl("positivity|stain", x, ignore.case = TRUE))), "YES", "NO")

progenes$Centrosome<-ifelse(progenes$Centrosome == "YES", 0.25, 0)
progenes$Cilia_expression<-ifelse(progenes$Cilia_expression == "YES", 0.5, 0)
progenes$Flagella<-ifelse(progenes$Flagella == "YES", 1, 0)
progenes$Cilia_localization<-ifelse(progenes$Cilia_localization == "YES", 1, 0)

progenes$score<-apply(progenes[,2:5], 1, max)
progenes<-hgncConverter(progenes, "Gene_name")  # This

protein_atlas_score<-progenes[,c(1,6)]

hgnc_names<-fread("hgnc_names.txt", select = c(1,3))
hgnc_names<-unique(hgnc_names)

protein_atlas_score<-left_join(hgnc_names[,2], protein_atlas_score, by = c("Approved symbol" = "Gene_name"))
protein_atlas_score$score[is.na(protein_atlas_score$score)]<-0
colnames(protein_atlas_score)[1:2]<-c("Gene_name", "pa_score")


write.table(protein_atlas_score, "proatlas_score_2022.txt", sep = "\t", row.names = FALSE, quote = FALSE)


procil<-semi_join(progenes, ciliaryGenes, by = "Gene_name")
noprocil<-semi_join(progenes, negativeCiliary, by = "Gene_name")


allprogenes<-left_join(unique(hgnc2[,1]), unique(progenes[,c(1,6)]), by = "Gene_name")

nontraining<-data.frame(Gene_name = training)
ciltraining<-ciliaryGenes[,1]
alltraining<-unique(rbind(ciltraining, nontraining))

alltraining$Gene_name[alltraining$Gene_name == "CILK"]<-"CILK1"
alltraining$Gene_name[alltraining$Gene_name == "GPBAR"]<-"GPBAR1"

allprogenes<-semi_join(allprogenes, alltraining, by = "Gene_name")


allprogenes$predicted<-ifelse(is.na(allprogenes$score), 0, 1)
allprogenes$is.ciliary<-ifelse(allprogenes$Gene_name %in% ciliaryGenes$Gene_name, 1, 0)
allprogenes$predicted<-as.factor(allprogenes$predicted)
allprogenes$is.ciliary<-as.factor(allprogenes$is.ciliary)
cnfs<-confusionMatrix(allprogenes$predicted, allprogenes$is.ciliary, positive = "1", mode="prec_recall")
cnfs$byClass[["F1"]]


proexp<-fread("expressionclustertissue_81_Ciliated.tsv")
colnames(proexp)[2]<-"Gene_name"

a<-semi_join(progenes, proexp, by = "Gene_name")

# Old codes

m1<-read_xlsx("modified_protein_atlas_gene_list.xlsx", sheet = 2)
m2<-read_xlsx("modified_protein_atlas_gene_list.xlsx", sheet = 3)
m3<-read_xlsx("modified_protein_atlas_gene_list.xlsx", sheet = 4)
m4<-read_xlsx("modified_protein_atlas_gene_list.xlsx", sheet = 5)

m1$score<-0.25
m2$score<-0.5
m3$score<-1
m4$score<-1

protein_atlas<-rbind(m1, m2, m3, m4)
protein_atlas<-subset(protein_atlas, !is.na(Gene_name))


write.table(protein_atlas[,c(1,12)], "protein_atlas_scores.csv", row.names = FALSE, quote = FALSE)




