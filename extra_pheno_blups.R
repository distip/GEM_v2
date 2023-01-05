library(Matrix)
library(lme4)
library(reshape2)
library(rrBLUP)
library(ggplot2)
library(tidyverse)
library(viridisLite)
library(viridis)
library(grid)

#################### Calculating blups for seperate N treatments ###########################3


spectra <- read.csv('BGEM_extra_phenotypes.csv')

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
str(spectra)

## removing block 4  ##
spectra[c(1051:1400), c('leaf_length', 'leaf_width' , 'ear_height', 'flag_leaf',  'plant_height') ] <- NA

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

bands <- colnames(spectra)[c(13:17)]
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

for(i in 2:2){
  for(j in 1:length(bands)){
    ### The hyperspectral data from just one N can be stored in a temporary variable.
    
    temp <- spectra.list[[i]]
    
    ### The hyperspectral data from just one band can be pulled out along with the field
    ### design variables. 
    
    temp <- temp[, which(colnames(temp) %in% c('genotype', 'Rep', 'Block', 'ASD', bands[j]))]
    colnames(temp)[5] <- 'reflectance'
    
    ### The BLUP model is 
    
    spectrum.blup.mod <- lmer(reflectance ~ (1|genotype) + (1|Rep) , data = temp)
    
    ### The variance components can be extracted to calculate broad-sense heritability.
    
    Vg <- data.frame(VarCorr(spectrum.blup.mod))$vcov[1]
    Ve <- data.frame(VarCorr(spectrum.blup.mod))$vcov[3]
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

ggplot(bands.H2[bands.H2$Trt == 'LN', ]) + 
  geom_bar(stat = 'identity' , aes(x= bands, y=H2))+
  labs(title = 'Broad Sense Heritabilities')+
  coord_cartesian(ylim=c(0.50, 1))+
  theme_bw(16)

names(spectra.blups.list[['HN']]) <- sub('.x', '', names(spectra.blups.list[['HN']]))

spectra_columns <- subset(spectra, select = c(1:12))
merged_1 <- merge(spectra_columns[which(spectra_columns$Trt== 'HN'),], spectra.blups.list[['HN']], by = 'genotype', all.x = TRUE)
merged_2 <- merge(spectra_columns[which(spectra_columns$Trt== 'LN'),], spectra.blups.list[['LN']], by = 'genotype', all.x = TRUE)


blups_merged <- merged_1 %>% full_join(merged_2)
blups_merged$Block <- factor(blups_merged$Block, levels= c('2', '4', '1', '3'))
write.csv(blups_merged, './extra_pheno_blups.csv', row.names = FALSE)
blups_merged <- read.csv('extra_pheno_blups.csv')

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


