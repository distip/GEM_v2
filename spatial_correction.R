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
spectra$rows <- as.numeric(spectra$rows)
spectra$ranges <- as.numeric(spectra$ranges)
spectra$PLOT.ID <- factor(spectra$PLOT.ID)

spectra <- subset(spectra, select = -c(X, Unnamed..0))


spectra.new <- spectra
spectra.new$Block <- factor(spectra.new$Block, levels= c('2', '4', '1', '3'))

ggplot(spectra.new, aes(rows, ranges, color=X730)) + 
  geom_point(size=1.3) +
  scale_y_continuous(name='ranges', limits = c(1,13))+
  geom_rect( aes(xmin=0.4, xmax = 51, ymin = 7.5, ymax = 12.5), fill=NA, colour='red')+
  #annotate('rect', xmin=0, xmax = 50, ymin = 7.5, ymax = 12.5, alpha= .1)+
  facet_wrap(Block~Trt)+
  scale_colour_viridis() +
  labs(title = 'Range and Row effects on leaf spectrum', caption = ' Blocks 1 and 4 = + N , 2 and 3 = -N\nred rectangles are the hybrids ')+
  theme_classic()
  #theme(strip.background = element_rect(color = 'black', fill = Trt))

spectra <- add_column(spectra, pos = numFactor(scale(as.numeric(spectra$rows)), scale(as.numeric(spectra$ranges))), .after = 'ranges' )
spectra <- add_column(spectra, ID = factor(rep(1, nrow(spectra))), .after= 'pos')

data <- spectra[which(spectra$Block == 4),]
data <- spectra
m1 <-glmmTMB(X730 ~ (1|genotype), data= data)
m2 <-glmmTMB(X730 ~ (1|genotype) + ar1(pos + 0 | ID) , data= data)
m3 <-glmmTMB(X730 ~ (1|genotype) + (1|ASD) + ar1(pos + 0 | ID), data=data)


anova(m1,m2,m3)



m4 <-glmmTMB(X750 ~ (1|Trt) + (1|genotype) + (1|ASD) + ar1(pos + 0 | ID), data= data)
m5 <-glmmTMB(X750 ~ (1|genotype) + (1|ASD) + (1|Trt/Block) + ar1(pos + 0 | ID), data= data)
spectra.nas <-spectra[which(is.na(spectra$X730)),]
