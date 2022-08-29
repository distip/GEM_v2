library(Matrix)
library(lme4)
library(reshape2)
library(rrBLUP)
library(ggplot2)
library(ggpubr)
library(tidyverse)
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


ggplot(data=calibration.melt[which(calibration.melt$ASD == '1'),], aes(band, reflectance)) + geom_line(color = 'green')

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

summary(anova)[[1]][['Pr(>F)']][2]

anova[2][1]

par(mfrow=c(2,2))
plot(anova)
par(mfrow=c(1,1))

two.way.plot <- ggplot(data=calibration.melt[which(calibration.melt$band == '900'),], aes(x=ASD, y=reflectance, groupr=PLOT.ID))+
  geom_point(cex=1.5, pch=1.0 , position=position_jitter(w=0.01, h=0))+
  geom_line(aes(group = PLOT.ID), size=0.2)
two.way.plot
