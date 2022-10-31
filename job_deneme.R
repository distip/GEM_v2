library(lme4)
library(reshape2)
library(rrBLUP)
library(ggplot2)
library(tidyverse)
library(zoom)
library(ggforce)
library(cowplot)


blups <- read.csv("spectra_blups.csv")

blups$Rep <- factor(blups$Rep)
blups$Block <- factor(blups$Block)
blups$year <- factor(blups$year)
blups$genotype <- factor(blups$genotype)
blups$note <- factor(blups$note)
blups$Trt <- factor(blups$Trt)
blups$ASD <- factor(blups$ASD)
blups$Group <- factor(blups$Group)
blups$rows <- factor(blups$rows)
blups$ranges <- factor(blups$ranges)
blups$PLOT.ID <- factor(blups$PLOT.ID)
blups$ASD  <- factor(blups$ASD)




########################   CALCULATION OF MID-PARENT HETEROSIS   #####################################
## heterosis calculation for low nitrogen condition
blups_LN <- blups[blups$Trt == "LN" & blups$Rep == 1 ,]

het <- blups_LN %>% filter(grepl(" X ", genotype))

mid_LN <- expand.grid(bands= 350:2500, genotype = het$genotype)
mid_LN$heterosis <- NA