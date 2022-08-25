
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

## Visualise for one genotype
spectra.sub <- spectra[which(spectra$genotype == 'BGEM-0292-S') ,]
spectra.sub.melt <- melt(spectra.sub, id.vars = c('PLOT.ID', 'Block' ,'Rep', 'Trt', 'year', 'genotype', 'note', 'Calibration', 'ASD'))

head(spectra.sub.melt)

colnames(spectra.sub.melt)[10:11] <- c('band', 'reflectance')

##  to plot the spectral data, the wavelength needs to be extracted from the band label by taking of the band number

spectra.sub.melt$band <- as.character(spectra.sub.melt$band)
spectra.sub.melt$band <- as.numeric(substr(spectra.sub.melt$band,2,5))

spectra.sub.melt$band

# since the field design has not been analyzed, there are three observations for this genotype
# we will plot only the first replicate.

ggplot(data=spectra.sub.melt, aes(band, reflectance, color=Trt)) + geom_line()

## to analyze the field design for spectral dataset, it is easiest to seperate the different treatments to create 
## an individual dataframe for each Trt. These can be stored in a list. Since there are 2 N app, the list length is 2 but we will
## put N in the model

## calculating the BLUPs for spectra
bands <- (colnames(spectra)[-c(1:9)])
bands_df <-as.data.frame(colnames(spectra)[-c(1:9)])
bands 
bands_df$H2 <- NA
colnames(bands_df) <- c('band', 'H2')
bands.H2 <- bands_df

spectral.blups.list <- data.frame(levels(spectra$genotype))
colnames(spectral.blups.list) <- 'genotype'



for(i in 1:length(bands)){
  temp <- spectra
  temp<-temp[, which(colnames(spectra) %in% c('genotype', 'Block', 'Trt', 'Rep', 'ASD', bands[i]))]
  colnames(temp)[6] <- 'reflectance'
  
  spectrum.blup.mod<-lmer(reflectance~(1|genotype)+(1|ASD)+(1|Trt:Rep)+(1|Block:Trt:Rep) , data=temp)

  Vg <- data.frame(VarCorr(spectrum.blup.mod))$vcov[1]
  Ve <- data.frame(VarCorr(spectrum.blup.mod))$vcov[5]
  H2 <- Vg/(Vg+(Ve/3))
  
  bands.H2[which(bands.H2$band==bands[i]), 'H2'] <- H2

  spectra.blups.temp <- ranef(spectrum.blup.mod)$genotype + fixef(spectrum.blup.mod)
  spectra.blups.temp <- data.frame(rownames(spectra.blups.temp), spectra.blups.temp)
  colnames(spectra.blups.temp) <- c('genotype', bands[i])
  
  spectral.blups.list <- merge(spectral.blups.list, spectra.blups.temp, by='genotype')
  print(i)
  
}


