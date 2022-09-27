
library(Matrix)
library(lme4)
library(reshape2)
library(rrBLUP)
library(ggplot2)
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

spectra <- subset(spectra, select = -c(X, Unnamed..0))

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
  
  spectrum.blup.mod<-lmer(reflectance~(1|genotype)+(1|Trt:Block), data=temp)

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

View(spectral.blups.list)
head(bands.H2)

bands.H2$band <- as.character(bands.H2$band)
bands.H2$band <- as.numeric(substr(bands.H2$band,2,5))

## plotting broad sense heritabilitiy

ggplot(bands.H2, aes(band ,H2)) +
  geom_line()+
  labs(title = 'Broad Sense Heritability of the Bands')

#################### Calculating blups for seperate N treatments ###########################3

levels(spectra$Trt)

spectra.list <- vector('list' , 2)
spectra.list

for(i in 1:2) {
  ##The appropriate N can be pulled out using which().
  
  spectra.list[[i]] <- spectra[which(spectra$Trt ==levels(spectra$Trt)[i]),]
  
  ### To keep track of where each N appears within the list, names() can be used to 
  ### set the name of each list element to the N
  
  names(spectra.list)[i] <- levels(spectra$Trt)[i]
  
}

View(spectra.list)

### A dataframe can be created to store the broad-sense heritability estimates for each
### band on each N. The expand.grid() function creates a data frame where each row 
### contains a different combination of the bands and N. An "H2" column can be added
### to store the broad-sense heritability estimates. 

bands <- colnames(spectra)[-c(1:11)]
bands

Trt <- as.character(levels(spectra$Trt))
Trt

bands.H2 <- expand.grid(bands, Trt)
bands.H2

colnames(bands.H2) <- c('bands', 'Trt')
bands.H2$H2 <- NA
head(bands.H2)

### Another list with a length of 2 can be created to store the BLUPs for the bands.

spectra.blups.list <- vector('list', 2)
spectra.blups.list

for(i in 1:2){
  
  ### To get started, a dataframe containing only the "gid" column can be created for each
  ### element of the list. 
  
  spectra.blups.list[[i]] <- data.frame(levels(spectra$genotype))
  colnames(spectra.blups.list[[i]]) <- 'genotype'
  
}
View(spectra.blups.list)

### As before, names() can be used to track the appropriate Ns.

names(spectra.blups.list) <- c('HN', 'LN')
### The following nested loop will first go through the list of dataframes each containing
### hyperspectral data from a single N.

for(i in 1:length(spectra.list)){
  for(j in 1:length(bands)){
    ### The hyperspectral data from just one N can be stored in a temporary variable.
    
    temp <- spectra.list[[i]]
    
    ### The hyperspectral data from just one band can be pulled out along with the field
    ### design variables. 
    
    temp <- temp[, which(colnames(temp) %in% c('genotype', 'Rep', 'Block', 'ASD', bands[j]))]
    colnames(temp)[5] <- 'reflectance'
    
    ### The BLUP model is 
    
    spectrum.blup.mod <- lmer(reflectance ~ (1|genotype)+(1|ASD)+(1|Block:Rep) , data = temp)
    
    ### The variance components can be extracted to calculate broad-sense heritability.
    
    Vg <- data.frame(VarCorr(spectrum.blup.mod))$vcov[1]
    Ve <- data.frame(VarCorr(spectrum.blup.mod))$vcov[4]
    H2 <- Vg/(Vg+(Ve/2))
    
    
    ### The broad-sense heritability can be put in the appropriate place within the 
    ### bands.H2 dataframe using which().
    
    bands.H2[which(bands.H2$bands == bands[j] & bands.H2$Trt == names(spectra.list)[i]), 'H2'] <- H2
    
    ### The BLUPs centered around the mean can also be calculated and stored in a dataframe. 
    
    spectra.blups.temp <- ranef(spectrum.blup.mod)$genotype + fixef(spectrum.blup.mod)
    spectra.blups.temp <- data.frame(rownames(spectra.blups.temp), spectra.blups.temp)
    colnames(spectra.blups.temp) <- c('genotype', bands[j])
    
    ### The BLUPs can be merged based on the "genotype" column with the existing dataframes 
    ### in the BLUPs list using the merge() function.
    
    spectra.blups.list[[i]] <- merge(spectra.blups.list[[i]], spectra.blups.temp, by = 'genotype')
    
    ### counters can be used to track the progress
    
    print(i)
    print(j)
    
  }
}

View(spectra)
head(bands.H2)


write.csv(bands.H2, './spectrum_heritabilities.csv', row.names = FALSE)
write.csv(spectra.blups.list[["HN"]], './spectrum_blups_HN.csv', row.names = FALSE)
write.csv(spectra.blups.list[["LN"]], './spectrum_blups_LN.csv', row.names = FALSE)


### To plot the broad-sense heritabilities for each band, the wavelengths need to be 
### extracted from the band labels as before.

bands.H2$bands <- as.character(bands.H2$bands)
bands.H2$bands<- as.numeric(substr(bands.H2$bands,2,5))

### To plot different lines for each of the different N application, the color can be set to "Trt".

ggplot(bands.H2, aes(bands, H2, color=Trt)) + geom_line()+
  labs(title = 'Broad Sense Heritabilities Under Different Nitrogen Applications')

names(spectra.blups.list[['HN']]) <- sub('.x', '', names(spectra.blups.list[['HN']]))

spectra_columns <- subset(spectra, select = c(1:11))
merged_1 <- merge(spectra_columns[which(spectra_columns$Trt== 'HN'),], spectra.blups.list[['HN']], by = 'genotype')
merged_2 <- merge(spectra_columns[which(spectra_columns$Trt== 'LN'),], spectra.blups.list[['LN']], by = 'genotype')


blups_merged <- merged_1 %>% full_join(merged_2)
blups_merged$Block <- factor(blups_merged$Block, levels= c('2', '4', '1', '3'))
write.csv(blups_merged, './spectra_blups.csv', row.names = FALSE)

ggplot(blups_merged, aes(rows, ranges, color=X730)) + 
  geom_point(size=1.3) +
  scale_y_continuous(name='ranges', limits = c(1,13))+
  geom_rect( aes(xmin=0.4, xmax = 51, ymin = 7.5, ymax = 12.5), fill=NA, colour='red')+
  #annotate('rect', xmin=0, xmax = 50, ymin = 7.5, ymax = 12.5, alpha= .1)+
  facet_wrap(~Trt)+
  scale_colour_viridis() +
  labs(title = 'Range and Row effects on leaf spectrum after spatial correction \n (1|genotype)+(1|ASD)+(1|Block:Rep)', caption = ' Blocks 1 and 4 = + N , 2 and 3 = -N\nred rectangles are the hybrids ')+
  theme_classic()
#theme(strip.background = element_rect(color = 'black', fill = Trt))


### To look at the hyperspectral data from multiple genotypes on a single data, the first 
### fifteen genotypes from the first date can be pulled out. As before, melt() moves the 
### spectral bands from columns to rows, and the wavelengths need to be extracted from the
### band labels.


spectra.blups.sub.melt <- melt(spectra.blups.list)
a <- sapply(strsplit(as.character(spectra.blups.sub.melt$variable), '[.]'), '[[', 1)
spectra.blups.sub.melt$wavelength <- as.numeric(substr(a,2,5))



View(spectra.blups.sub.melt)
View(spectra.blups.list)
colnames(spectra.blups.sub.melt)[4] <- 'Trt'


write.csv(spectra.blups.sub.melt, './fullspectra_plots.csv', row.names = FALSE )
spectra.blups.sub.melt <- read.csv('fullspectra_plots.csv')

ggplot(spectra.blups.sub.melt[which(spectra.blups.sub.melt$Trt == 'HN'),], aes(wavelength, value, color=genotype)) +
  geom_line(size=0.3)+
  theme(legend.position = 'none')


ggplot() +
  geom_line(data=spectra.blups.sub.melt[which(spectra.blups.sub.melt$Trt == 'LN'),],aes(wavelength, value, group= genotype, color='LN'), size=0.1)+
  geom_line(data=spectra.blups.sub.melt[which(spectra.blups.sub.melt$Trt == 'HN'),], aes(wavelength, value, group= genotype, color='HN'), size=0.1)+
  scale_color_manual(name='Treatment',  values = c('HN'= 'green', 'LN' = 'yellow'))
  

ggplot(data = spectra.blups.sub.melt)+
  geom_line(aes(wavelength, value, group=genotype, color=Trt), size=0.1)+
  labs(title = 'Leaf spectra (BLUPS)', x='bands', y='reflectance')



######### Calculating blups for different ASDs #########
spectra <- read.csv("Raw_spectrum_merged")
spectra <- spectra[which(spectra$ASD %in% c(1,2)),]
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

levels(spectra$ASD)

spectra.list <- vector('list' , 2)
spectra.list

for(i in 1:2) {
  ##The appropriate asd can be pulled out using which().
  
  spectra.list[[i]] <- spectra[which(spectra$ASD ==levels(spectra$ASD)[i]),]
  
  ### To keep track of where each N appears within the list, names() can be used to 
  ### set the name of each list element to the N
  
  names(spectra.list)[i] <- levels(spectra$ASD)[i]
  
}

View(spectra.list)

### A dataframe can be created to store the broad-sense heritability estimates for each
### band on each N. The expand.grid() function creates a data frame where each row 
### contains a different combination of the bands and N. An "H2" column can be added
### to store the broad-sense heritability estimates. 

bands <- colnames(spectra)[-c(1:9)]
bands

ASD <- as.character(levels(spectra$ASD))
ASD

bands.H2 <- expand.grid(bands, ASD)
bands.H2

colnames(bands.H2) <- c('bands', 'ASD')
bands.H2$H2 <- NA
head(bands.H2)

### Another list with a length of 2 can be created to store the BLUPs for the bands.

spectra.blups.list <- vector('list', 2)
spectra.blups.list

for(i in 1:2){
  
  ### To get started, a dataframe containing only the "gid" column can be created for each
  ### element of the list. 
  
  spectra.blups.list[[i]] <- data.frame(levels(spectra$genotype))
  colnames(spectra.blups.list[[i]]) <- 'genotype'
  
}
View(spectra.blups.list)

### As before, names() can be used to track the appropriate Ns.

names(spectra.blups.list) <- c('ASD1', 'ASD2')
### The following nested loop will first go through the list of dataframes each containing
### hyperspectral data from a single N.

for(i in 1:length(spectra.list)){
  for(j in 1:length(bands)){
    ### The hyperspectral data from just one N can be stored in a temporary variable.
    
    temp <- spectra.list[[i]]
    
    ### The hyperspectral data from just one band can be pulled out along with the field
    ### design variables. 
    
    temp <- temp[, which(colnames(temp) %in% c('genotype', 'Rep', 'Block', 'Trt', bands[j]))]
    colnames(temp)[5] <- 'reflectance'
    
    ### The BLUP model is 
    
    spectrum.blup.mod <- lmer(reflectance ~ (1|genotype)+(1|Trt)+(1|Rep:Trt)+(1|Rep:Trt:Block), data = temp)
    
    ### The variance components can be extracted to calculate broad-sense heritability.
    
    Vg <- data.frame(VarCorr(spectrum.blup.mod))$vcov[1]
    Ve <- data.frame(VarCorr(spectrum.blup.mod))$vcov[5]
    H2 <- Vg/(Vg+(Ve/2))
    
    
    ### The broad-sense heritability can be put in the appropriate place within the 
    ### bands.H2 dataframe using which().
    
    bands.H2[which(bands.H2$bands == bands[j] & bands.H2$ASD == names(spectra.list)[i]), 'H2'] <- H2
    
    ### The BLUPs centered around the mean can also be calculated and stored in a dataframe. 
    
    spectra.blups.temp <- ranef(spectrum.blup.mod)$genotype + fixef(spectrum.blup.mod)
    spectra.blups.temp <- data.frame(rownames(spectra.blups.temp), spectra.blups.temp)
    colnames(spectra.blups.temp) <- c('genotype', bands[j])
    
    ### The BLUPs can be merged based on the "genotype" column with the existing dataframes 
    ### in the BLUPs list using the merge() function.
    
    spectra.blups.list[[i]] <- merge(spectra.blups.list[[i]], spectra.blups.temp, by = 'genotype')
    
    ### counters can be used to track the progress
    
    print(i)
    print(j)
    
  }
}

View(spectra)
head(bands.H2)


write.csv(bands.H2, './spectrum_heritabilities_forASDsS.csv', row.names = FALSE)
write.csv(spectra.blups.list[["HN"]], './spectrum_blups_HN.csv', row.names = FALSE)
write.csv(spectra.blups.list[["LN"]], './spectrum_blups_LN.csv', row.names = FALSE)


### To plot the broad-sense heritabilities for each band, the wavelengths need to be 
### extracted from the band labels as before.

bands.H2$bands <- as.character(bands.H2$bands)
bands.H2$bands<- as.numeric(substr(bands.H2$bands,2,5))

### To plot different lines for each of the different N application, the color can be set to "Trt".

ggplot(bands.H2, aes(bands, H2, color=ASD)) + geom_line()+
  labs(title = 'Broad Sense Heritabilities Under Different Nitrogen Applications')
