
#################theKaplanLab#############

#################Hasan Can Demirci############


############### Necessary Libraries ###########

library(ggplot2)
library(dplyr)

## working directory ##

setwd("C:\\Users\\Hasan Can Demirci\\Desktop\\CilioGenics Revision\\histograms\\hst_table5s")


####################### Uploading the File ####################

hcd <- read.csv2("Table_S5_Comparative.xlsx - table_S3.csv", sep = ",", header = TRUE)

hcd <- unique(hcd)


## summarize data to create difference between type ##

hcd_processed <- hcd %>%
  group_by(Ciliary, Single_cell) %>%
  summarize(count = n()) %>%
  ungroup()

## get plot ##

ggplot(hcd_processed, aes(x = Ciliary, y = count, fill = Single_cell)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_text(aes(label = count), position = position_dodge(width = 0.9), vjust = -0.25) +
  scale_fill_manual(values = c("YES" = "#DDCC77", "NO" = "#CC6677")) +
  ylim(0, 110) +
  ggtitle("Cluster 31 and Cluster 37") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(),
        plot.title = element_text(hjust = 0.5))

### saving plot ###

ggsave("s_table5s.tiff", dpi = 300, width = 11, height = 7)

