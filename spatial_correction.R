library(Matrix)
library(lme4)
library(reshape2)
library(rrBLUP)
library(ggplot2)
library(glmmTMB)
library(ncf)
library(viridisLite)
library(viridis)
library(tibble)
library(PerformanceAnalytics)

spectra <- read.csv("Raw_spectrum_merged")
View(spectra)

str(spectra)

spectra$Rep <- factor(spectra$Rep)
spectra$Block <- factor(spectra$Block)
spectra$year <- factor(spectra$year)
spectra$genotype <- factor(spectra$genotype)
spectra$note <- factor(spectra$note)
spectra$Trt <- factor(spectra$Trt)
spectra$ASD <- factor(spectra$ASD)
spectra$rows <- factor(spectra$rows)
spectra$ranges <- factor(spectra$ranges)
spectra$PLOT.ID <- factor(spectra$PLOT.ID)

spectra <- subset(spectra, select = -c(X, Unnamed..0))

spectra <- add_column(spectra, pos = numFactor(scale(as.numeric(spectra$rows)), scale(as.numeric(spectra$ranges))), .after = 'ranges' )
spectra <- add_column(spectra, ID = factor(rep(1, nrow(spectra))), .after= 'pos')

data <- spectra[which(spectra$Block == 3),]
data <- spectra
m1 <-glmmTMB(X1400 ~ (1|genotype), data= data)
m2 <-glmmTMB(X1400 ~ (1|genotype) + (1|ASD) , data= data)
m3 <-glmmTMB(X1400 ~ (1|genotype) + (1|ASD) + ar1(pos + 0 | ID), data=data)
m4 <-glmmTMB(X1400 ~ (1|Trt) + (1|genotype) + (1|ASD) + ar1(pos + 0 | ID), data= data)
m5 <-glmmTMB(X1400 ~ (1|genotype) + (1|ASD) + (1|Trt/Block) + ar1(pos + 0 | ID), data= data)

anova(m1,m2,m3)
