library(Matrix)
library(lme4)
library(reshape2)
library(rrBLUP)
library(ggplot2)
library(tidyverse)
library(viridisLite)
library(viridis)
library(grid)



data <- read_csv('extra_pheno_blues.csv')

View(data)

melt = melt(data=data, id =c("genotype"   ,  "PLOT.ID"    ,  "rows"   ,      "ranges"   ,    "Block"    ,    "Rep"        ,  "Group"  ,      "Trt"    ,      "year" ,        "note" ,       
                             "Calibration",  "ASD" ))
View(melt)


ggplot(data=melt, aes(x=Trt , y=value)) + 
  geom_boxplot()+
  facet_wrap(vars(variable) , scales = 'free')+
  theme_bw(16)
