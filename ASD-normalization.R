spectra <- read.csv('Raw_spectrum_merged')

str(spectra)
spectra <- spectra[which(spectra$ASD %in% c('1','2')),]

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

levels(spectra$ASD)
ASD.list <- vector('list', 2)
ASD.list

for(i in 1:2){
  ASD.list[[i]] <- spectra[which(spectra$ASD == levels(spectra$ASD)[i]),]
  names(ASD.list)[i] <- levels(spectra$ASD)[i]
}
View(ASD.list)


spectra.stat.list <- vector('list', 2)

for (i in 1:2){
  spectra.stat.list[[i]] <- data.frame(bands= 350:2500)
}

names(spectra.stat.list) <- c(1,2)
View(spectra.stat.list)

df <- data.frame(ASD= '', median = '', mean = '', sd = '')

for(j in 1:2 ){
  for(i in 12:length(colnames(spectra))){
    IiM <- median(spectra[which(spectra$ASD == j ), i])
    Iim <- mean(spectra[which(spectra$ASD == j), i])
    Iisd <- sd(spectra[which(spectra$ASD == j ), i])
    #df <- data.frame(median = IiM, mean = Iim, sd = Iisd)
    df <- rbind(df, data.frame(ASD = j, median = IiM, mean = Iim, sd = Iisd))
    #spectra.stat.list[[j]] <- merge(spectra.stat.list[[j]], df, by= 'bands' )  
  }
}
df.stats <- df[-1,]
df.stats$bands <- rep(350:2500, 2) 

df.stats$median <- as.numeric(df.stats$median)
df.stats$mean <- as.numeric(df.stats$mean)
df.stats$sd <- as.numeric(df.stats$sd)


for(j in 1:2){
  for (i in 12:length(colnames(spectra))){
    normalized <- (spectra[which(spectra$ASD==j), i] - df.stats[which(df.stats$ASD == j),][(i-11), 'mean'])/df.stats[which(df.stats$ASD == j),][(i-11), 'sd']
    ASD.list[[j]][,i] <- normalized
  }
}

View(ASD.list)

scaled <- do.call('rbind', ASD.list)


#################### Calculating blups for seperate N treatments ###########################3
library(Matrix)
library(lme4)
library(reshape2)
library(rrBLUP)
library(ggplot2)



spectra <- scaled
str(spectra)

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
    
    spectrum.blup.mod <- lmer(reflectance ~ (1|genotype)+(1|Block)+(1|Block:Rep), data = temp)
    
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
  labs(title = 'Broad Sense Heritabilities Under Different Nitrogen Applications (after scale-median)')




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

