
library(Matrix)
library(lme4)
library(reshape2)
library(rrBLUP)
library(ggplot2)
library(tidyverse)
library(viridisLite)
library(viridis)



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
spectra$Group <- factor(spectra$Group)
spectra$rows <- factor(spectra$rows)
spectra$ranges <- factor(spectra$ranges)
spectra$PLOT.ID <- factor(spectra$PLOT.ID)
spectra$ASD  <- factor(spectra$ASD)
spectra$Calibration <- factor(spectra$Calibration)

spectra <- subset(spectra, select = -c(X))

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
    
    spectrum.blup.mod <- lmer(reflectance ~ (1|genotype)+(1|ASD)+(1|Rep) , data = temp)
    
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
  labs(title = 'Broad Sense Heritabilities Under Different Nitrogen Applications')+
  theme_bw()

names(spectra.blups.list[['HN']]) <- sub('.x', '', names(spectra.blups.list[['HN']]))

spectra_columns <- subset(spectra, select = c(1:11))
merged_1 <- merge(spectra_columns[which(spectra_columns$Trt== 'HN'),], spectra.blups.list[['HN']], by = 'genotype')
merged_2 <- merge(spectra_columns[which(spectra_columns$Trt== 'LN'),], spectra.blups.list[['LN']], by = 'genotype')


blups_merged <- merged_1 %>% full_join(merged_2)
blups_merged$Block <- factor(blups_merged$Block, levels= c('2', '4', '1', '3'))
write.csv(blups_merged, './spectra_blups.csv', row.names = FALSE)
blups_merged <- read.csv('spectra_blups.csv')

str(blups_merged)

blups_merged$genotype <- factor(blups_merged$genotype)
blups_merged$PLOT.ID <- factor(blups_merged$PLOT.ID)
blups_merged$rows <- factor(blups_merged$rows)
blups_merged$ranges <- factor(blups_merged$ranges)
blups_merged$Block <- factor(blups_merged$Block)
blups_merged$Rep <- factor(blups_merged$Rep)
blups_merged$Trt  <- factor(blups_merged$Trt)
blups_merged$year <- factor(blups_merged$year)
blups_merged$note <- factor(blups_merged$note)
blups_merged$Calibration <- factor(blups_merged$Calibration)
blups_merged$ASD <- factor(blups_merged$ASD)

blups_merged$note <- NA

for(i in 1:length(blups_merged$genotype)){
  gen <- as.character(blups_merged$genotype[i])
  if(grepl(' X ', gen)){
    blups_merged$note[i] <- 'Hybrid'
  }
  else{
    blups_merged$note[i] <- 'Inbred'
  }
}

spectra.blups.melt.2 <- melt(blups_merged, id.vars = c('PLOT.ID', 'genotype','rows', 'ranges', 'Block' ,'Rep', 'Trt', 'year', 'note', 'Calibration', 'ASD'))
a <- sapply(strsplit(as.character(spectra.blups.melt.2$variable), '[.]'), '[[', 1)
spectra.blups.melt.2$wavelength <- as.numeric(substr(a,2,5))

data <- spectra.blups.melt.2[spectra.blups.melt.2$Trt == 'HN',] %>% group_by(wavelength,note) %>% 
  summarise(mean.ref = mean(value, na.rm=TRUE), sd.ref = sd(value, na.rm = TRUE), se.ref= sd(value, na.rm=TRUE)/sqrt(length(value)), 
            max = max(value, na.rm = TRUE), min = min(value, na.rm = TRUE))


ggplot(blups_merged, aes(rows, ranges, color=X730)) + 
  geom_point(size=1.3) +
  scale_y_continuous(name='ranges', limits = c(1,13))+
  geom_rect( aes(xmin=0.4, xmax = 51, ymin = 7.5, ymax = 12.5), fill=NA, colour='red')+
  #annotate('rect', xmin=0, xmax = 50, ymin = 7.5, ymax = 12.5, alpha= .1)+
  facet_wrap(~Trt)+
  scale_colour_viridis() +
  labs(title = 'Spatial Map of leaf spectrum after BLUPs \n (1|genotype)+(1|ASD)+(1|Rep)', caption = ' Blocks 1 and 4 = + N , 2 and 3 = -N\nred rectangles represent the hybrids ')+
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

data <- spectra.blups.sub.melt %>% group_by(wavelength,Trt) %>% 
  summarise(mean.ref = mean(value, na.rm=TRUE), sd.ref = sd(value, na.rm = TRUE), se.ref= sd(value, na.rm=TRUE)/sqrt(length(value)), 
            max = max(value, na.rm = TRUE), min = min(value, na.rm = TRUE))

data$wavelength <- as.numeric(data$wavelength)

plt_blups <- ggplot(data=data, aes(x=wavelength, group= note)) +
  geom_line(aes(y=mean.ref, color=note), size = 0.6)+
  geom_ribbon(aes(ymin=mean.ref-se.ref , ymax=mean.ref+se.ref , fill=note),alpha=0.3)+
  labs(title = 'Leaf Spectra inbred vs hybrids under HN', caption = '**Envelopes represent 1 sd from the mean')+
  theme_bw()

plt_blups



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





                                 ######################    CALCULATING BLUPS FOR FOR EACH UNIQUE COMBINATIONS ######################

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
spectra$Group <- factor(spectra$Group)
spectra$rows <- factor(spectra$rows)
spectra$ranges <- factor(spectra$ranges)
spectra$PLOT.ID <- factor(spectra$PLOT.ID)
spectra$ASD  <- factor(spectra$ASD)
spectra$Calibration <- factor(spectra$Calibration)

spectra <- subset(spectra, select = -c(Unnamed..0, X))

new_GID <- paste(spectra$genotype, spectra$Trt, sep= "_")

spectra_comb <- add_column(spectra, new_GID = new_GID, .after='PLOT.ID')

str(spectra_comb)

spectra_comb$Rep <- factor(spectra_comb$Rep)
spectra_comb$Block <- factor(spectra_comb$Block)
spectra_comb$year <- factor(spectra_comb$year)
spectra_comb$genotype <- factor(spectra_comb$genotype)
spectra_comb$note <- factor(spectra_comb$note)
spectra_comb$Trt <- factor(spectra_comb$Trt)
spectra_comb$ASD <- factor(spectra_comb$ASD)
spectra_comb$Group <- factor(spectra_comb$Group)
spectra_comb$rows <- factor(spectra_comb$rows)
spectra_comb$ranges <- factor(spectra_comb$ranges)
spectra_comb$PLOT.ID <- factor(spectra_comb$PLOT.ID)
spectra_comb$ASD  <- factor(spectra_comb$ASD)
spectra_comb$Calibration <- factor(spectra_comb$Calibration)
spectra_comb$new_GID  <- factor(spectra_comb$new_GID)
levels(spectra_comb$Trt)


levels(spectra_comb$Trt)

spectra_comb.list<- spectra_comb

View(spectra_comb.list)

### A dataframe can be created to store the broad-sense heritability estimates for each
### band on each N. The expand.grid() function creates a data frame where each row 
### contains a different combination of the bands and N. An "H2" column can be added
### to store the broad-sense heritability estimates. 

bands <- colnames(spectra_comb)[-c(1:13)]
bands

Trt <- as.character(levels(spectra_comb$Trt))
Trt

bands.H2 <- expand.grid(bands)
bands.H2

colnames(bands.H2) <- c('bands')
bands.H2$H2 <- NA
head(bands.H2)

### Another list with a length of 2 can be created to store the BLUPs for the bands.

spectra_comb.blups.list <- data.frame(levels(spectra_comb$new_GID))
colnames(spectra_comb.blups.list) <- c('new_GID')
View(spectra_comb.blups.list)


### The following nested loop will first go through the list of dataframes each containing
### hyperspectra_combl data from a single N.


for(j in 1:length(bands)){
### The hyperspectra_combl data from just one N can be stored in a temporary variable.
    
  temp <- spectra_comb.list
    
  ### The hyperspectra_combl data from just one band can be pulled out along with the field
  ### design variables. 
    
  temp <- temp[, which(colnames(temp) %in% c('new_GID', 'Rep', 'Block', 'ASD', bands[j]))]
  colnames(temp)[5] <- 'reflectance'
    
  ### The BLUP model is 
    
  spectrum.blup.mod <- lmer(reflectance ~ (1|new_GID) + (1|ASD) + (1|Rep) , data = temp)
    
  ### The variance components can be extracted to calculate broad-sense heritability.
    
  Vg <- data.frame(VarCorr(spectrum.blup.mod))$vcov[1]
  Ve <- data.frame(VarCorr(spectrum.blup.mod))$vcov[4]
  H2 <- Vg/(Vg+(Ve/2))
    
    
  ### The broad-sense heritability can be put in the appropriate place within the 
  ### bands.H2 dataframe using which().
    
  bands.H2[which(bands.H2$bands == bands[j]), 'H2'] <- H2
    
  ### The BLUPs centered around the mean can also be calculated and stored in a dataframe. 
    
  spectra_comb.blups.temp <- ranef(spectrum.blup.mod)$new_GID + fixef(spectrum.blup.mod)
  spectra_comb.blups.temp <- data.frame(rownames(spectra_comb.blups.temp), spectra_comb.blups.temp)
  colnames(spectra_comb.blups.temp) <- c('new_GID', bands[j])
    
  ### The BLUPs can be merged based on the "genotype" column with the existing dataframes 
  ### in the BLUPs list using the merge() function.
    
  spectra_comb.blups.list <- merge(spectra_comb.blups.list, spectra_comb.blups.temp, by = 'new_GID')
    
  ### counters can be used to track the progress
    
  print(j)
    
}

### To plot the broad-sense heritabilities for each band, the wavelengths need to be 
### extracted from the band labels as before.

bands.H2$bands <- as.character(bands.H2$bands)
bands.H2$bands<- as.numeric(substr(bands.H2$bands,2,5))

### To plot different lines for each of the different N application, the color can be set to "Trt".

ggplot(bands.H2, aes(bands, H2)) + geom_line()+
  labs(title = 'Broad Sense Heritabilities Under Different Nitrogen Applications')+
  theme_bw()


names(spectra_comb.blups.list) <- sub('.x', '', names(spectra_comb.blups.list))

spectra_comb_columns <- subset(spectra_comb, select = c(1:13))
merged_1 <- merge(spectra_comb_columns, spectra_comb.blups.list, by = 'new_GID')


blups_merged_v2 <- merged_1
blups_merged_v2$Block <- factor(blups_merged$Block, levels= c('2', '4', '1', '3'))
write.csv(blups_merged_v2, './spectra_comb_blups_v2.csv', row.names = FALSE)
blups_merged_v2 <- read.csv('spectra_comb_blups_v2.csv')

str(blups_merged_v2)

blups_merged_v2$genotype <- factor(blups_merged_v2$genotype)
blups_merged_v2$PLOT.ID <- factor(blups_merged_v2$PLOT.ID)
blups_merged_v2$rows <- factor(blups_merged_v2$rows)
blups_merged_v2$ranges <- factor(blups_merged_v2$ranges)
blups_merged_v2$Block <- factor(blups_merged_v2$Block)
blups_merged_v2$Rep <- factor(blups_merged_v2$Rep)
blups_merged_v2$Trt  <- factor(blups_merged_v2$Trt)
blups_merged_v2$year <- factor(blups_merged_v2$year)
blups_merged_v2$note <- factor(blups_merged_v2$note)
blups_merged_v2$Calibration <- factor(blups_merged_v2$Calibration)
blups_merged_v2$ASD <- factor(blups_merged_v2$ASD)
blups_merged_v2$new_GID <- factor(blups_merged_v2$new_GID)

blups_merged_v2 <- blups_merged_v2[blups_merged_v2$Rep == 2 , ]


data <- spectra_comb.blups.melt.2[spectra_comb.blups.melt.2$Trt == 'HN',] %>% group_by(wavelength,note) %>% 
  summarise(mean.ref = mean(value, na.rm=TRUE), sd.ref = sd(value, na.rm = TRUE), se.ref= sd(value, na.rm=TRUE)/sqrt(length(value)), 
            max = max(value, na.rm = TRUE), min = min(value, na.rm = TRUE))


ggplot(blups_merged_v2, aes(rows, ranges, color=X550)) + 
  geom_point(size=1.3) +
  scale_y_continuous(name='ranges', limits = c(1,13))+
  geom_rect( aes(xmin=0.4, xmax = 51, ymin = 7.5, ymax = 12.5), fill=NA, colour='red')+
  #annotate('rect', xmin=0, xmax = 50, ymin = 7.5, ymax = 12.5, alpha= .1)+
  facet_wrap(~Trt)+
  scale_colour_viridis() +
  labs(title = 'Spatial Map of leaf spectrum after BLUPs \n (1|genotype)+(1|ASD)+(1|Rep)', caption = ' Blocks 1 and 4 = + N , 2 and 3 = -N\nred rectangles represent the hybrids ')+
  theme_classic()
#theme(strip.background = element_rect(color = 'black', fill = Trt))


### To look at the hyperspectra_combl data from multiple genotypes on a single data, the first 
### fifteen genotypes from the first date can be pulled out. As before, melt() moves the 
### spectra_combl bands from columns to rows, and the wavelengths need to be extracted from the
### band labels.


spectra_comb.blups.sub.melt <- melt(spectra_comb.blups.list)
a <- sapply(strsplit(as.character(spectra_comb.blups.sub.melt$variable), '[.]'), '[[', 1)
spectra_comb.blups.sub.melt$wavelength <- as.numeric(substr(a,2,5))



View(spectra_comb.blups.sub.melt)
View(spectra_comb.blups.list)
colnames(spectra_comb.blups.sub.melt)[4] <- 'Trt'


write.csv(spectra_comb.blups.sub.melt, './fullspectra_comb_plots.csv', row.names = FALSE )
spectra_comb.blups.sub.melt <- read.csv('fullspectra_comb_plots.csv')

data2 <- data %>% group_by(wavelength,Trt) %>% 
  summarise(mean.ref = mean(value, na.rm=TRUE), sd.ref = sd(value, na.rm = TRUE), se.ref= sd(value, na.rm=TRUE)/sqrt(length(value)), 
            max = max(value, na.rm = TRUE), min = min(value, na.rm = TRUE))


data <-melt(blups_merged_v2, id.vars = c('PLOT.ID','new_GID', 'genotype','rows', 'ranges', 'Block' ,'Rep', 'Trt', 'year', 'note', 'Group', 'Calibration', 'ASD'))
a <- sapply(strsplit(as.character(data$variable), '[.]'), '[[', 1)
data$wavelength <- as.numeric(substr(a,2,5))
View(data)


data$wavelength <- as.numeric(data$wavelength)

plt_blups <- ggplot(data=data2, aes(x=wavelength, group= Trt)) +
  geom_line(aes(y=mean.ref, color=Trt), size = 0.6)+
  geom_ribbon(aes(ymin=mean.ref-se.ref , ymax=mean.ref+se.ref),alpha=0.3)+
  labs(title = 'Leaf spectra_comb inbred vs hybrids under HN', caption = '**Envelopes represent 1 sd from the mean')+
  theme_bw()

plt_blups



ggplot(spectra_comb.blups.sub.melt[which(spectra_comb.blups.sub.melt$Trt == 'HN'),], aes(wavelength, value, color=genotype)) +
  geom_line(size=0.3)+
  theme(legend.position = 'none')


ggplot() +
  geom_line(data=data[which(data$Trt == 'LN'),],aes(wavelength, value, group= new_GID, color='LN'), size=0.1)+
  geom_line(data=data[which(data$Trt == 'HN'),], aes(wavelength, value, group= new_GID, color='HN'), size=0.1)+
  scale_color_manual(name='Treatment',  values = c('HN'= 'green', 'LN' = 'yellow'))


ggplot(data = data)+
  geom_line(aes(wavelength, value, group=new_GID, color=Trt), size=0.1)+
  labs(title = 'Leaf spectra_comb (BLUPS)', x='bands', y='reflectance')+
  theme_bw()




