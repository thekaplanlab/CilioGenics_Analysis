#################theKaplanLab#############

#################ComplexUpset & Observation of Intersection Members ############

############ Hasan Can Demirci #######
              
############### Necessary Libraries ###########

library(ggplot2)
library(ComplexUpset)
library(readxl)
library(dplyr)
library(tibble)
library(writexl)

####################### Uploading the File ####################

#setwd("/Users/oktayismailkaplan/Desktop/Manuscripts/Ciliogenics")
hcd <- read_xlsx("/Users/oktayismailkaplan/Desktop/Manuscripts/Ciliogenics/Figures/TF_genes/TF_genes1.xlsx")
Genes = as.data.frame(hcd)

####################### Collecting of All Genes under a column #################

Genes <- c(hcd$RFX2, hcd$RFX3, hcd$MYB, hcd$GLIS3, hcd$JAZF1, hcd$FOXJ)
Genes
####################### Making new Table including Genes and T. Factors ########

Genes_table <- tibble(Genes)
Genes_table
Genes_table <- data.frame(Genes_table, RFX2= 0, RFX3= 0,
                          MYB= 0,GLIS3= 0,JAZF1= 0,FOXJ= 0)
###################### Determine the values based on existence #################

Genes_table$RFX2 <- ifelse(Genes_table$Genes %in% hcd$RFX2, 1, 0)
Genes_table$RFX3 <- ifelse(Genes_table$Genes %in% hcd$RFX3, 1, 0)
Genes_table$MYB <- ifelse(Genes_table$Genes %in% hcd$MYB, 1, 0)
Genes_table$GLIS3 <- ifelse(Genes_table$Genes %in% hcd$GLIS3, 1, 0)
Genes_table$JAZF1 <- ifelse(Genes_table$Genes %in% hcd$JAZF1, 1, 0)
Genes_table$FOXJ <- ifelse(Genes_table$Genes %in% hcd$FOXJ, 1, 0)
Genes_table
Genes_table <- na.omit(Genes_table)

##################### Determination of the intersection group name #############

tf = colnames(Genes_table)[2:7]
tf
##################### Determination of Set_size function to arrange set size ###

set_size = function(w, h, factor=1.5) {
  s = 1 * factor
  options(
    repr.plot.width=w * s,
    repr.plot.height=h * s,
    repr.plot.res=100 / factor,
    jupyter.plot_mimetypes='image/png',
    jupyter.plot_scale=1
  )
}

set_size(2386,6)

#######################Determination of the query_by_degree to give properties##
####################### Determination of Intersection properties################

query_by_degree = function(data, groups, params_by_degree, shared, ...) {
  intersections = unique(upset_data(data, groups, ...)$plot_intersections_subset)
  lapply(
    intersections,
    FUN=function(x) {
      members = ComplexUpset:::get_intersection_members(x)[[1]]
      degree = as.character(ComplexUpset:::calculate_degree(x))
      
      args = c(
        list(intersect=members),
        shared,
        params_by_degree[[degree]]
      )
      do.call(upset_query, args)
    }
  )
}

######################### Calling ComplexUpset Plot ############################
upset(
  Genes_table,
  tf,
  width_ratio=0.1,
  min_size = 15,
  sort_intersections_by=c('degree', 'cardinality'),
  matrix=intersection_matrix(
    geom=geom_point(shape='circle filled', size=3)
  ),
      queries=
        query_by_degree(
          Genes_table,
          tf,
          min_size = 15,
          params_by_degree=list(
            '1'=list(fill='#CCFFFF'),
            '2'=list(fill='#CCE5FF'),
            '3'=list(fill='#CCCCFF'),
            '4'=list(fill='#E5CCFF'),
            '5'=list(fill='#FFCCFF'),
            '6'=list(fill='#FFCCE5')
          ),
          shared=list(
            only_components=c("intersections_matrix", "Intersection size"),
            color='black'
          )
        )
)

############################### Making Distinctsssss ###########################

distinct_genes_table <- distinct(Genes_table,Genes,RFX2, RFX3, MYB,GLIS3,JAZF1,FOXJ)

upset(
  distinct_genes_table,
  tf,
  width_ratio=0.1,
  min_size = 15,
  sort_intersections_by=c('degree', 'cardinality'),
  matrix=intersection_matrix(
    geom=geom_point(shape='circle filled', size=3)
  ),
  queries=
    query_by_degree(
      distinct_genes_table,
      tf,
      min_size = 15,
      params_by_degree=list(
        '1'=list(fill='#CCFFFF'),
        '2'=list(fill='#CCE5FF'),
        '3'=list(fill='#CCCCFF'),
        '4'=list(fill='#E5CCFF'),
        '5'=list(fill='#FFCCFF'),
        '6'=list(fill='#FFCCE5')
      ),
      shared=list(
        only_components=c("intersections_matrix", "Intersection size"),
        color='black'
      )
    )
)


upset(
  distinct_genes_table,
  tf,
  width_ratio=0.1,
  min_size = 15,
  sort_intersections_by=c('degree', 'cardinality'),
  matrix=intersection_matrix(
    geom=geom_point(shape='circle filled', size=3)
  ),
  queries=
    query_by_degree(
      distinct_genes_table,
      tf,
      min_size = 15,
      params_by_degree=list(
        '1'=list(fill='#6D5271'),
        '2'=list(fill='#FFC77D'),
        '3'=list(fill='#FCE7CC'),
        '4'=list(fill='#DCD3E2'),
        '5'=list(fill='#616D3D'),
        '6'=list(fill='#429EBD')
      ),
      shared=list(
        only_components=c("intersections_matrix", "Intersection size"),
        color='black'
      )
    )
)

#################### Finding Intersection Members from Upset_Data ##############

upset_data <- upset_data(distinct_genes_table, tf)
upset_data
################# Finding Distincsss again #####################################

distinct1 <- distinct(upset_data$presence,Genes, intersection)

######################## Making Group and count the members number##############

last <- distinct1 %>%
  group_by(distinct1$intersection) %>%
  summarize(Num_Genes = n(), Genes = paste(Genes, collapse = ", "))

############################ Creation of Excel File on Destkop #################

write_xlsx(last, "Cmplx_Upset_Members.xlsx")
