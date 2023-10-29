#################theKaplanLab#############

#################Hasan Can Demirci############

#################ComplexUpset & Observation of Intersection Members ############



############### Necessary Libraries ###########

library(ggplot2)
library(ComplexUpset)
library(readxl)
library(dplyr)
library(tibble)
library(writexl)
library(stringr)
####################### Uploading the File ####################

setwd("C:\\Users\\Hasan Can Demirci\\Desktop\\CilioGenics Revision\\histograms")

hcd <- read.csv2("all_genes_scores.csv")

## select required columns ##

hcd <- hcd %>%
  select(1,3,4,5,6,8)

hcd <- unique(hcd)

## create table only genes ##

all_genes <- hcd[1]
all_genes <- unique(all_genes)

colnames(hcd)

## Create the data for complexUpset -First columns creation ###

last_data <- data.frame(all_genes, Reyfman = "",
                        Carraro = "", Habermann = "", Murthy = "",
                        Trachea = "")

hcd[is.na(hcd)] <- ""
###################### Determine the values based on existence #################

last_data$Reyfman <- ifelse(hcd$Reyfman >=0, 1, 0)
last_data$Carraro <- ifelse(hcd$Carraro>=0, 1, 0)
last_data$Habermann <- ifelse(hcd$Habermann>=0, 1, 0)
last_data$Murthy <- ifelse(hcd$Murthy>=0, 1, 0)
last_data$Trachea <- ifelse( hcd$Trachea>=0, 1, 0)


last_data <- na.omit(last_data)

##################### Determination of the intersection group name #############

tf = colnames(last_data)[2:6]

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

distinct_last_data <- distinct(last_data,Gene_name,Reyfman, Carraro, Habermann,Murthy,Trachea)

## remove all columns == 0 ##

distinct_last_data <- distinct_last_data[rowSums(distinct_last_data[, 2:6]) != 0, ]

set_size(4756,5)
upset(
  distinct_last_data,
  tf,
  set_sizes=FALSE,
  width_ratio=0.1,
  min_size = 15,
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
      distinct_last_data,
      tf,
      min_size = 15,
      params_by_degree=list( ### colors for color-blind-friendly ###
        '1'=list(fill='#632D5A'),
        '2'=list(fill='#882255'),
        '3'=list(fill='#AA4499'),
        '4'=list(fill='#CC6677'),
        '5'=list(fill='#DDCC77'),
        '6'=list(fill='#117733'),
        '7'=list(fill='#44AA99')
      ),
      shared=list(
        only_components=c("intersections_matrix", "Intersection size"),
        color='black'
      )
    )
)

### saving plot ###

ggsave("s_table2s.tiff", dpi = 300, width = 11, height = 7)

