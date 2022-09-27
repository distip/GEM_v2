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

spectra <- subset(spectra, select = -c(X, Unnamed..0))


#################### Calculating blues for seperate N treatments ###########################3

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


### Another list with a length of 2 can be created to store the BLUPs for the bands.

spectra.blues.list <- vector('list', 2)
spectra.blues.list

for(i in 1:2){
  
  ### To get started, a dataframe containing only the "gid" column can be created for each
  ### element of the list. 
  
  spectra.blues.list[[i]] <- data.frame(levels(spectra$genotype))
  colnames(spectra.blues.list[[i]]) <- 'genotype'
  
}
View(spectra.blues.list)

### As before, names() can be used to track the appropriate Ns.

names(spectra.blues.list) <- c('HN', 'LN')
### The following nested loop will first go through the list of dataframes each containing
### hyperspectral data from a single N.

for(i in 2:length(spectra.list)){
  for(j in 1:length(bands)){
    ### The hyperspectral data from just one N can be stored in a temporary variable.
    
    temp <- spectra.list[[i]]
    
    ### The hyperspectral data from just one band can be pulled out along with the field
    ### design variables. 
    
    temp <- temp[, which(colnames(temp) %in% c('genotype', 'Rep', 'Block', 'ASD', bands[j]))]
    colnames(temp)[5] <- 'reflectance'
    
    ### The BLUE model is 
    
  spectrum.blue.mod <- lmer(reflectance ~ genotype + (1|ASD) + (1|Rep) , REML = FALSE,  data = temp)
    
  
     ### The BLUEs centered around the mean can also be calculated and stored in a dataframe. 
    spectra.int <- fixef(spectrum.blue.mod)[1]
    spectra.blues <- fixef(spectrum.blue.mod)
    #spectra.blues[-1] <- spectra.blues[-1] + spectra.int
    spectra.blues.temp <- spectra.blues[-1] + spectra.int
    spectra.blues.temp <- data.frame (spectra.blues.temp)
    spectra.blues.temp<- cbind(genotype=rownames(spectra.blues.temp), spectra.blues.temp)
    rownames(spectra.blues.temp) <- NULL
    spectra.blues.temp$genotype <- gsub('genotype', '', spectra.blues.temp$genotype)
    colnames(spectra.blues.temp) <- c('genotype', bands[j])
    
    ### The BLUPs can be merged based on the "genotype" column with the existing dataframes 
    ### in the BLUPs list using the merge() function.
    
    spectra.blues.list[[i]] <- merge(spectra.blues.list[[i]], spectra.blues.temp, by = 'genotype')
    
    ### counters can be used to track the progress
    
    print(i)
    print(j)
    
  }
}

View(spectra)
head(bands.H2)



### To plot the broad-sense heritabilities for each band, the wavelengths need to be 
### extracted from the band labels as before.

bands.H2$bands <- as.character(bands.H2$bands)
bands.H2$bands<- as.numeric(substr(bands.H2$bands,2,5))

### To plot different lines for each of the different N application, the color can be set to "Trt".

ggplot(bands.H2, aes(bands, H2, color=Trt)) + geom_line()+
  labs(title = 'Broad Sense Heritabilities Under Different Nitrogen Applications')

names(spectra.blups.list[['HN']]) <- sub('.x', '', names(spectra.blups.list[['HN']]))

spectra_columns <- subset(spectra, select = c(1:11))
merged_1 <- merge(spectra_columns[which(spectra_columns$Trt== 'HN'),], spectra.blues.list[['HN']], by = 'genotype')
merged_2 <- merge(spectra_columns[which(spectra_columns$Trt== 'LN'),], spectra.blues.list[['LN']], by = 'genotype')


blues_merged <- merged_1 %>% full_join(merged_2)
blues_merged$Block <- factor(blues_merged$Block, levels= c('2', '4', '1', '3'))
write.csv(blues_merged, './spectra_blues.csv', row.names = FALSE)

ggplot(blues_merged, aes(rows, ranges, color=X730)) + 
  geom_point(size=1.3) +
  scale_y_continuous(name='ranges', limits = c(1,13))+
  geom_rect( aes(xmin=0.4, xmax = 51, ymin = 7.5, ymax = 12.5), fill=NA, colour='red')+
  #annotate('rect', xmin=0, xmax = 50, ymin = 7.5, ymax = 12.5, alpha= .1)+
  facet_wrap(~Trt)+
  scale_colour_viridis() +
  labs(title = 'BLUEs Corrected \n genotype+(1|ASD)+(1|Block:Rep)', caption = ' Blocks 1 and 4 = + N , 2 and 3 = -N\nred rectangles are the hybrids ')+
  theme_classic()
#theme(strip.background = element_rect(color = 'black', fill = Trt))


### To look at the hyperspectral data from multiple genotypes on a single data, the first 
### fifteen genotypes from the first date can be pulled out. As before, melt() moves the 
### spectral bands from columns to rows, and the wavelengths need to be extracted from the
### band labels.


spectra.blues.sub.melt <- melt(spectra.blues.list)
a <- sapply(strsplit(as.character(spectra.blues.sub.melt$variable), '[.]'), '[[', 1)
spectra.blues.sub.melt$wavelength <- as.numeric(substr(a,2,5))



View(spectra.blues.sub.melt)
View(spectra.blues.list)
colnames(spectra.blues.sub.melt)[4] <- 'Trt'


write.csv(spectra.blues.sub.melt, './fullspectra_plots_blues.csv', row.names = FALSE )
spectra.blues.sub.melt <- read.csv('fullspectra_plots_blues.csv')

ggplot(spectra.blues.sub.melt[which(spectra.blues.sub.melt$Trt == 'HN'),], aes(wavelength, value, color=genotype)) +
  geom_line(size=0.3)+
  theme(legend.position = 'none')


ggplot() +
  geom_line(data=spectra.blues.sub.melt[which(spectra.blues.sub.melt$Trt == 'LN'),],aes(wavelength, value, group= genotype, color='LN'), size=0.1)+
  geom_line(data=spectra.blues.sub.melt[which(spectra.blues.sub.melt$Trt == 'HN'),], aes(wavelength, value, group= genotype, color='HN'), size=0.1)+
  scale_color_manual(name='Treatment',  values = c('HN'= 'green', 'LN' = 'yellow'))


ggplot(data = spectra.blues.sub.melt)+
  geom_line(aes(wavelength, value, group=genotype, color=Trt), size=0.1)+
  labs(title = 'Leaf spectra (BLUPS)', x='bands', y='reflectance')
