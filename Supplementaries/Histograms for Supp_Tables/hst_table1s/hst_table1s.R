
#################theKaplanLab#############

#################Hasan Can Demirci############


############### Necessary Libraries ###########

library(ggplot2)

## working directory ##

setwd("C:\\Users\\Hasan Can Demirci\\Desktop\\CilioGenics Revision\\histograms\\hst_table1s")


####################### Uploading the File ####################

hcd <- read.csv2("Table_S1.xlsx - Table_S1A.csv",sep = ",")

hcd <- unique(hcd)

## add count colon with counted number ##

hcd$count <- ifelse(hcd$Ciliary== "YES", 91, 503)

## to ne frame to add text on the bars ##

hcd <- data.frame(
  Ciliary = c("YES", "Unknown"),
  count = c(91, 503)
)


## get plot ##

ggplot(hcd, aes(x = Ciliary, y= count)) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = Ciliary)) +
  geom_text(aes(label = count), position = position_dodge(width=0.9), vjust = -0.25) +
  scale_fill_manual(values = c("YES" = "#CC6677", "Unknown" = "#DDCC77")) +
  ylim(0, 505) +
  ggtitle("Cao et. al. (2017)") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(),
        plot.title = element_text(hjust = 0.5))

### saving plot ###

ggsave("s_table1s.tiff", dpi = 300, width = 11, height = 7)

