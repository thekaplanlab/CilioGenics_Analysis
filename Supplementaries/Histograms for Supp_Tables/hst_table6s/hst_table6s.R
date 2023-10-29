
#################theKaplanLab#############

#################Hasan Can Demirci############


############### Necessary Libraries ###########

library(ggplot2)
library(dplyr)

## working directory ##

setwd("C:\\Users\\Hasan Can Demirci\\Desktop\\CilioGenics Revision\\histograms\\hst_table6s")


####################### Uploading the File ####################

hcd <- read.csv2("Table_S6_TF.xlsx - Table_S4.csv",sep = ",")

hcd <- unique(hcd)

## summurazi data ## 

hcd_processed <- hcd %>%
  group_by(Ciliary.Gene) %>%
  summarize(count = n()) %>%
  ungroup() 

## get plot ##

ggplot(hcd_processed, aes(x = Ciliary.Gene, y= count)) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = Ciliary.Gene)) +
  geom_text(aes(label = count), position = position_dodge(width=0.9), vjust = -0.25) +
  scale_fill_manual(values = c("YES" = "#CC6677", "Unknown" = "#DDCC77")) +
  ylim(0, 1305) +
  ggtitle("") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(),
        plot.title = element_text(hjust = 0.5))

### saving plot ###

ggsave("s_table6s.tiff", dpi = 300, width = 11, height = 7)

