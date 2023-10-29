#################theKaplanLab#############

#################Hasan Can Demirci############



############### Necessary Libraries ###########

library(ggplot2)
library(ComplexUpset)
library(readxl)
library(dplyr)
library(tibble)
library(writexl)
library(stringr)
####################### Uploading the File ####################

setwd("C:\\Users\\Hasan Can Demirci\\Desktop\\CilioGenics Revision\\histograms\\hst_table8s")

hcd <- read.csv2("Table_S8_Protein_atlas.xlsx - Table_S6.csv",sep = ",", header = T)

## select required columns ##

hcd <- hcd %>%
  select(1,2,3,4,5)

hcd <- unique(hcd)

###################### Determine the values based on existence #################

hcd$Centrosome <- ifelse(hcd$Centrosome =="YES", 1, 0)
hcd$Cilia_expression <- ifelse(hcd$Cilia_expression =="YES", 1, 0)
hcd$Flagella <- ifelse(hcd$Flagella =="YES", 1, 0)
hcd$Cilia_localization <- ifelse(hcd$Cilia_localization =="YES", 1, 0)


hcd <- na.omit(hcd)

##################### Determination of the intersection group name #############

tf = colnames(hcd)[2:5]

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

### Make distinct to reach single correct data ###

distinct_hcd <- distinct(hcd,Gene_name, Centrosome ,Cilia_expression, Flagella, Cilia_localization)

## remove all columns == 0 ##

distinct_hcd <- distinct_hcd[rowSums(distinct_hcd[, 2:5]) != 0, ]

set_size(584,4)

## do not forget to arrange min_size to see all intersection ##

upset(
  distinct_hcd,
  tf,
  set_sizes=FALSE,
  width_ratio=0.1,
  min_size = 2,
  sort_intersections_by=c('degree', 'cardinality'),
  
  themes=upset_modify_themes( ### arrangement of intersection table with dots ###
    
    list('intersections_matrix'=theme(axis.text.y=element_text(size=20), axis.text.x = element_text(size = 0),
                                      axis.title.x.bottom =element_text(size=20) , axis.text=element_text(size=20, face='bold'),
                                      axis.title=element_text(size=20, face='bold'),
                                      axis.title.y =element_text(size=50), axis.title.y.left = element_text(size=5),
                                      axis.text.y.left =element_text(size=20), axis.title.x.top = element_text(size=20),
                                      axis.ticks.y = element_line(size=1), axis.line.y.left = element_line(size=1),
                                      axis.line.x =element_line(size=1) ),
         
         'Intersection size'=theme( ### arrangement of intersection table with bars ###
           
           axis.text=element_text(size=20, face='bold'),
           axis.title=element_text(size=20)
         )
    )
  ),
  matrix=intersection_matrix(
    geom=geom_point(shape='circle filled', size=3),
  ),
  queries=
    query_by_degree(
      distinct_hcd,
      tf,
      min_size = 2,
      params_by_degree=list( ### colors for color-blind-friendly ###
        '1'=list(fill='#632D5A'),
        '2'=list(fill='#882255'),
        '3'=list(fill='#AA4499'),
        '4'=list(fill='#CC6677'),
        '5'=list(fill='#DDCC77')
      ),
      shared=list(
        only_components=c("intersections_matrix", "Intersection size"),
        color='black'
      )
    )
)
### saving plot ###

ggsave("s_table8s.tiff", dpi = 300, width = 11, height = 7)

