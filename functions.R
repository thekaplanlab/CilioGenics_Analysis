



`%fin%` <- function(x, table) {
  stopifnot(require(fastmatch))
  fmatch(x, table, nomatch = 0L) > 0L
}

hgncConverter2<-function (genelist, colname) 
{
  colname2 <- "Gene_name"
  colname3 <- "Genename"
  tempcolname <- "Gene_synonyms"
  approved <- semi_join(genelist, hgnc, by = setNames(colname2, 
                                                      colname))
  someWOapproved <- anti_join(genelist, approved, by = colname)
  notApproved <- semi_join(someWOapproved, hgnc, by = setNames(tempcolname, 
                                                               colname))
  some <- rbind(approved, notApproved)
  non <- anti_join(genelist, some, by = colname)
  if ("Gene_name" %in% colnames(genelist) & !("Genename" %in% 
                                              colnames(genelist))) {
    hgnc_temp <- hgnc
    colnames(hgnc_temp)[1] <- "Genename"
    temp_notApproved <- left_join(notApproved, hgnc_temp, 
                                  by = setNames(tempcolname, colname))
    temp_notApproved[[colname]] <- temp_notApproved$Genename
    temp_notApproved <- temp_notApproved %>% dplyr::select(!Genename)
  }
  else if (!("Gene_name" %in% colnames(genelist))) {
    temp_notApproved <- left_join(notApproved, hgnc, by = setNames(tempcolname, 
                                                                   colname))
    temp_notApproved[[colname]] <- temp_notApproved$Gene_name
    temp_notApproved <- temp_notApproved %>% dplyr::select(!Gene_name)
  }
  else {
    stop("Change the name of that damn column!")
  }
  all <- rbind(approved, temp_notApproved, non)
  return(all)
}


mousesynonymConverter<-function(genelist,colname){
  
  names(mouse_synonyms1)[1]<-"Gene_name_temp"
  n_occur <- data.frame(table(mouse_synonyms1$Gene_synonyms))
  bg5<-mouse_synonyms1[mouse_synonyms1$Gene_synonyms %in% n_occur$Var1[n_occur$Freq > 1],]
  
  bg5<-bg5[which(bg5$Gene_synonyms %in% genelist[[colname]][which(!(genelist[[colname]] %in% mouse_synonyms1$Gene_name))]),]
  colnames(bg5)[2]<-colname
  
  genelist<-merge(genelist, bg5, by = colname, all = TRUE, allow.cartesian = TRUE)
  for (i in 1:length(genelist[[colname]])){
    
    if (!is.na(genelist$Gene_name_temp[i])){
      genelist[[colname]][i]<-genelist$Gene_name_temp[i]
    }
  }
  pb = txtProgressBar(min = 0, max =length(genelist[[colname]]) , initial = 0)
  for (i in 1:length(genelist[[colname]])){
    if (!(genelist[[colname]][i] %fin% mouse_synonyms1$Gene_name_temp) && length(mouse_synonyms1$Gene_name_temp[which(mouse_synonyms1$Gene_synonyms %fin% genelist[[colname]][i])]) == 1){
      genelist[[colname]][i]<-mouse_synonyms1$Gene_name_temp[which(mouse_synonyms1$Gene_synonyms %fin% genelist[[colname]][i])]
    }
    setTxtProgressBar(pb,i)
  }
  return(genelist[,1:(length(genelist)-1), drop = FALSE])
}


mousegnameConverter<-function(genelist,colname){
  
  names(mouse_homology)[1]<-"Gene_name_temp"
  n_occur <- data.frame(table(mouse_homology$Gene_synonyms))
  bg5<-mouse_homology[mouse_homology$Gene_synonyms %in% n_occur$Var1[n_occur$Freq > 1],]
  
  bg5<-bg5[which(bg5$Gene_synonyms %in% genelist[[colname]][which(!(genelist[[colname]] %in% mouse_homology$Gene_name))]),]
  colnames(bg5)[2]<-"Gene_name_A"
  
  genelist<-merge(genelist, bg5, by = colname, all = TRUE, allow.cartesian = TRUE)
  for (i in 1:length(genelist[[colname]])){
    
    if (!is.na(genelist$Gene_name_temp[i])){
      genelist[[colname]][i]<-genelist$Gene_name_temp[i]
    }
  }
  
  pb = txtProgressBar(min = 0, max =length(genelist[[colname]]) , initial = 0)
  for (i in 1:length(genelist[[colname]])){
    if (!(genelist[[colname]][i] %fin% mouse_homology$Gene_name_temp) && length(mouse_homology$Gene_name_temp[which(mouse_homology$Gene_synonyms %fin% genelist[[colname]][i])]) == 1){
      genelist[[colname]][i]<-mouse_homology$Gene_name_temp[which(mouse_homology$Gene_synonyms %fin% genelist[[colname]][i])]
    }
    setTxtProgressBar(pb,i)
  }
  return(genelist[,1:(length(genelist)-1), drop = FALSE])
}



apiToGene<- function(x){
  dt<-data.frame(t(x), stringsAsFactors = FALSE)
  dt<-as.data.frame(dt[grepl("symbol", dt[[1]]),])
  dt<-as.data.frame(str_split(dt[[1]], pattern = '"symbol":"', simplify = TRUE))
  dt$V2<-gsub('"', "", dt$V2)
  dt<-as.data.frame(dt[,2])
  colnames(dt)<-"gene_name"
  dt
}

normalization0<-function(x){
  a<-(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))
}

normalization<-function(x){
  if(length(x[!is.na(x)]) > 0){
    a<-(x-0)/(max(x, na.rm = TRUE)-0)
  }
  else {a<-NA}
}

normalization2<-function(x){
  if(length(x[!is.na(x)]) > 0){
    a<-(x-0)/(max(x, na.rm = TRUE)-0)
  }
  else {a<-NA}
}


stdize = function(x, ...) {((x - min(x, ...)) / (max(x, ...) - min(x, ...)))-1}

outliers<-function(x){
  qnt <- quantile(x, probs=c(.25, .75), na.rm = T)
  caps <- quantile(x, probs=c(.05, .95), na.rm = T)
  H <- 3 * IQR(x, na.rm = T)
  ifelse(x > (qnt[2] + H) & !is.na(x), caps[2], x)
  # x[which(x < (qnt[1] - H) & !is.na(x))] <- caps[1]
  # x[which(x > (qnt[2] + H) & !is.na(x))] <- caps[2]
}

impute.mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))
scoring<-function(data){
  data %>% filter(`Experimental System Type` == "physical") %>%
    mutate(Score = gsub("-", NA, Score)) %>%
    mutate(Score = as.numeric(Score)) %>%
    group_by(`Publication Source`) %>%
    mutate(ppi_score = outliers(Score)) %>%
    mutate(ppi_score = normalization(ppi_score)) %>%
    ungroup() %>%
    mutate(ppi_score = impute.mean(ppi_score)) %>%
    distinct()
}

filtering<-function(file_location, e, bitscore){
  
  downloads<-read_xlsx("files/sequence_downloads.xlsx")
  downloads$link<-str_extract(downloads$link, "[^\\/]+(?=_protein\\.faa.gz$)")
  
  # Read all files
  files<- list.files(file_location)
  all<-data.frame(human_np = NA)
  
  for (i in 1:length(files)){
    
    # Name the columns of all files and merge them ----
    aa<-fread(paste0(file_location, files[i]), select = c(1,2,11,12)) %>%
      arrange(desc(V12)) %>%
      distinct(V1, .keep_all = TRUE)
    aa<-aa %>% filter(V11 < e & V12 > bitscore) %>% dplyr::select(1,2)
    nm<-gsub(".txt","",files[i])
    spname<-downloads$species[downloads$link == nm]
    colnames(aa)[1]<-"human_np"
    colnames(aa)[2]<-spname
    all<-merge(all, aa, by = "human_np", all = TRUE)
  }
  
  # Score genes based on presence or absence in organisms ----
  all1<-as.data.frame(sapply(all[,2:length(all)], function(x) ifelse(is.na(x), 0, 1)))
  all1<-cbind(all[,1],all1)
  colnames(all1)[1]<-"human_np"
  
  # Change order of species in merged table ----
  orderoforg<-read.table("files/organism_order.txt", sep = "\t", header = TRUE)
  orgorder<-intersect(orderoforg$Organism, colnames(all1))
  
  all1<-all1[,c(colnames(all1)[1],orgorder)]
  all1[,2:length(all1)] <- lapply(all1[,2:length(all1)], factor)
  
  # Annotate protein ids with gene names ----
  annot<-fread("files/human_annot.txt") %>% dplyr::select(Locus, `Protein product`) %>%
    unique()
  colnames(annot)<-c("Gene_name","human_np")
  annot$human_np<-str_split(annot$human_np, "\\.", simplify = TRUE)[,1]
  
  all1<-left_join(all1, annot, by = "human_np")
  all1<-all1[!duplicated(all1$Gene_name),]
  
  newscores<-all1
  colnames(newscores)<-gsub("Â","",colnames(newscores))
  colnames(newscores)<-gsub(" ","_",colnames(newscores))
  colnames(newscores)<-gsub(strsplit(colnames(newscores)[3],"")[[1]][4],"_",colnames(newscores))
  newscores<-newscores[,2:74]
  newscores[,1:72] <- lapply(newscores[,1:72], factor)
  newscores<-hgncConverter(newscores, "Gene_name")
  newscores<-newscores[-1,]
  newscores<-unique(newscores)
  
  return(newscores)
}

