library(Matrix)
library(lme4)
library(reshape2)
library(rrBLUP)
library(ggplot2)
library(tidyverse)
library(viridisLite)
library(viridis)
library(grid)

data <- read_csv('blups_merged_v2')
data <- data.frame(data)

spectra$Rep <- factor(spectra$Rep)
spectra$Block <- factor(spectra$Block)
spectra$year <- factor(spectra$year)
spectra$genotype <- factor(spectra$genotype)
spectra$note <- factor(spectra$note)
spectra$Trt <- factor(spectra$Trt)
spectra$ASD <- factor(spectra$ASD)
spectra$Group <- factor(spectra$Group)
spectra$rows <- factor(spectra$rows)
spectra$ranges <- factor(spectra$ranges)
spectra$PLOT.ID <- factor(spectra$PLOT.ID)
spectra$ASD  <- factor(spectra$ASD)
spectra$Calibration <- factor(spectra$Calibration)

View(data)

x <- data[data['Trt'] == 'HN', 'X550'] 
y <- data[data['Trt'] == 'LN', 'X550']

fit.lm <- lm(x ~ y)

slope <- coef(fit.lm)[2]


ggplot(data= data %>% select(genotype,X550,Trt))+
  geom_point(aes(x=Trt, y=X550))+
  geom_line(aes(x=Trt, y=X550, group=genotype, color='red'))

             