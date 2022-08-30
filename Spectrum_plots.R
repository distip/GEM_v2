library(Matrix)
library(lme4)
library(reshape2)
library(rrBLUP)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(MASS)
library(heritability)
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
spectra <- subset(spectra, select = -X)

spectra.sub <- spectra #[which(spectra$ASD %in% c(1,2)),]
spectra.sub.melt <- melt(spectra.sub, id.vars = c('PLOT.ID', 'Block' ,'Rep', 'Trt', 'year', 'genotype', 'note', 'Calibration', 'ASD'))
View(spectra.sub.melt)

colnames(spectra.sub.melt)[10:11] <- c('band', 'reflectance')

spectra.sub.melt$band <- as.character(spectra.sub.melt$band)
spectra.sub.melt$band <- as.numeric(substr(spectra.sub.melt$band,2,5))

spectra.sub.melt$band

str(spectra.sub.melt)

spectra.sub.melt.grouped <- spectra.sub.melt %>% group_by(genotype, band, Trt) %>% summarise(reflect=mean(reflectance, na.rm = TRUE))

str(spectra.sub.melt.grouped)

plot.son <- ggplot() + 
geom_line(data=spectra.sub.melt.grouped[which(spectra.sub.melt.grouped$Trt== 'HN'),], aes(band, reflect, group=genotype,  color='HN'), size=0.2)+ 
  geom_line(data=spectra.sub.melt.grouped[which(spectra.sub.melt.grouped$Trt== 'LN'),], aes(band, reflect, group=genotype,  color='LN'), size=0.2, alpha=0.4) +
  scale_colour_manual(name= 'Treatment', values=c('HN' = 'green', 'LN' = 'yellow'))+
  labs(title = 'Leaf spectra under two nitrogen treatments (raw spectrum)')
  

plot.son

ggplot(data=spectra.sub.melt.grouped) + geom_line(aes(band, reflect, color=Trt, alpha=0.4))+
  labs(title = 'Leaf Spectra (Raw data)', x='bands', y='reflectance')



##### Plotting the calibrations

calibration <- read.csv('calibration')
View(calibration)

colnames(calibration)
calibration <- subset(calibration, select = -c(X))

calibration.melt <- melt(calibration, id.vars = c('PLOT.ID', 'ASD'))
View(calibration.melt)

colnames(calibration.melt)[3:4] <- c('band', 'reflectance')

calibration.melt$band <- as.character(calibration.melt$band)
calibration.melt$band <- as.numeric(substr(calibration.melt$band, 2, 5))

calibration.melt$band

calibration.melt$ASD <-factor(calibration.melt$ASD)

str(calibration.melt)

calibration.melt$PLOT.ID <- factor(calibration.melt$PLOT.ID)
calibration.melt$band <- factor(calibration.melt$band)


ggplot() + 
  geom_line(data=calibration.melt[which(calibration.melt$ASD == '2'),],aes(band,reflectance, group=PLOT.ID, color='ASD1'),size=0.2)+
  geom_line(data=calibration.melt[which(calibration.melt$ASD == '1'),],aes(band,reflectance, group=PLOT.ID, color='ASD2'),size=0.2)+
  scale_colour_manual(name= 'ASD', values=c('ASD1' = 'green', 'ASD2' = 'yellow'))+
  labs(title = 'calibration results')


ggplot(data=calibration.melt[which(calibration.melt$ASD == '1'),], aes(band, reflectance)) + geom_line(color = 'green')


##### ANOVA ######


pvals <- c()
for(i in 350:2500) {
  anova <- aov(reflectance ~ PLOT.ID * ASD , data = calibration.melt[which(calibration.melt$band == i),])
  p <- summary(anova)[[1]][['Pr(>F)']][2]
  pvals <- c(pvals, p)
}

pvals.df <- data.frame(c(350:2500), pvals)
colnames(pvals.df)[1] <- 'band'


ggplot(pvals.df) + geom_point(aes(band, pvals), colour='black', size=1, shape=1) + geom_hline(yintercept = 0.05, linetype='dashed', color='red') +
  labs(title = 'p-values from ANOVA - comparing two ASDs')+ 
  scale_y_continuous(breaks= c(0.00, 0.05, 0.25, 0.50, 0.75, 1.00))


anova <- aov(reflectance ~ PLOT.ID * ASD , data = calibration.melt[which(calibration.melt$band == '2200'),])
summary(anova)

summary(anova)[[1]][['Mean Sq']][4]

anova[2][1]

par(mfrow=c(2,2))
plot(anova)
par(mfrow=c(1,1))

two.way.plot <- ggplot(data=calibration.melt[which(calibration.melt$band %in% c('350', '500', '750', '1000', '1250', '1500', '1750', '2000','2250','2500')),], aes(x=ASD, y=reflectance, group=PLOT.ID))+
  geom_point(cex=1.5, pch=1.0 , position=position_jitter(w=0.01, h=0))+
  geom_line(aes(group = PLOT.ID), size=0.2)+
  facet_grid (. ~ band)+
  labs(title = 'Leaf Reflectance obtained from 2 instruments in 10 bands', x='ASDs')
two.way.plot




##### Heritability ######

H_all <- c()
for(i in 3:ncol(calibration)) {
  H <- repeatability(data.vector= calibration[,i], geno.vector = calibration$PLOT.ID, covariates.frame = as.data.frame(calibration[,2]))
  H_all <- c(H_all, H$repeatability)
}

H2_df <- data.frame(c(350:2500), H_all)
colnames(H2_df)[1] <- 'band'

ggplot(data = H2_df) + 
  geom_line(aes(band, H_all)) + 
  labs(title = 'Repeatability of two instruments with the same leaf samples (40 leaves)')

