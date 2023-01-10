library(tidyverse)
library(ggplot2)
library(Matrix)
library(reshape2)



data <- read.csv('Spectrum_with_biochemical_traits_outlier_removal.csv')
colnames(data)[0:25]

data_sub <- subset(data, select = c('PLOT.ID', "X..N" ,"X..P" ,"X..K", "X..S", "X..Ca", "X..Mg", "ppm.Zn", "ppm.Fe", "ppm.Mn",  "ppm.Cu", 
                                    "ppm.B", "ppm.Mo", "Rep", "Group", "Trt", "genotype" ,"note"))
View(data_sub)

data_non_na <- data_sub
View(data_non_na)

melted <- melt(data_non_na, id.vars= c('PLOT.ID', 'Group', 'Trt', 'genotype', 'note', 'Rep'))
View(melted)


ggplot(data = melted[melted$Group =='Inbred',])  +
  geom_density(aes(value, fill=Trt), alpha=0.3)+
  theme_bw(14)+
  facet_wrap(vars(variable), scales = 'free')+
  theme(axis.ticks= element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())+
  labs(title= 'Inbreds under HN and LN - After outlier (1.5 IQR)')
