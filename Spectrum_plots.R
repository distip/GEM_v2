library(Matrix)
library(lme4)
library(reshape2)
library(rrBLUP)
library(ggplot2)

spectra <- read.csv("Raw_spectrum_merged")
View(spectra)

str(spectra)

spectra$Rep <- factor(spectra$Rep)
spectra$Block <- factor(spectra$Block)
spectra$year <- factor(spectra$year)
spectra$genotype <- factor(spectra$genotype)
spectra$note <- factor(spectra$note)
spectra$Trt <- factor(spectra$Trt)

spectra <- subset(spectra, select = -X)

spectra.sub <- spectra[which(spectra$ASD %in% c(1,2)),]
spectra.sub.melt <- melt(spectra.sub, id.vars = c('PLOT.ID', 'Block' ,'Rep', 'Trt', 'year', 'genotype', 'note', 'Calibration', 'ASD'))
View(spectra.sub.melt)

colnames(spectra.sub.melt)[10:11] <- c('band', 'reflectance')

spectra.sub.melt$band <- as.character(spectra.sub.melt$band)
spectra.sub.melt$band <- as.numeric(substr(spectra.sub.melt$band,2,5))

spectra.sub.melt$band

ggplot(data=spectra.sub.melt[which(spectra.sub.melt$Trt == 'LN'),], aes(band, reflectance, color=ASD)) + geom_line()



##### Plotting the calibrations

calibration <- read.csv('calibration')
View(calibration)

colnames(calibration)
calibration <- subset(calibration, select = -c(X, Unnamed..0 ))

calibration.melt <- melt(calibration, id.vars = c('PLOT.ID', 'ASD'))
View(calibration.melt)

colnames(calibration.melt)[3:4] <- c('band', 'reflectance')

calibration.melt$band <- as.character(calibration.melt$band)
calibration.melt$band <- as.numeric(substr(calibration.melt$band, 2, 5))

calibration.melt$band


ggplot(data=calibration.melt[which(calibration.melt$ASD == 2),], aes(band, reflectance, color= ASD)) + geom_line()


