library(tibble)
library(tidyr)
library(writexl)
library(dplyr)


a <- read.table("/Users/oktayismailkaplan/Desktop/Manuscripts/Ciliogenics/Figures/Publications/publications.txt",
                header = FALSE, sep = "|")

a
p <- a$V2 #Assign column name
g <- a$V1 #Assign column name



gp <- cbind(g,p) #combine two of them

gp <- gp[-1,] #remove the numbers in the first column
g
gp <- as.data.frame(gp)  #organize data into a 2-dimensional table of rows and columns,
gp

#summarize is a part of library(dplyr) function

gp_summary <- gp %>% #calculate the number of genes among publications
  group_by(g) %>%
  summarize(Num_Genes = n(), 
            p_values = toString(unique(p)))
gp_summary 
write_xlsx(gp_summary, "gp_summary.xlsx")



