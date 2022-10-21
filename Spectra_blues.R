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


### Another list with a length of 2 can be created to store the BLUEs for the bands.

spectra.blues.list <- vector('list', 2)
spectra.blues.list

for(i in 1:2){
  
  ### To get started, a dataframe containing only the "gid" column can be created for each
  ### element of the list. 
  
  spectra.blues.list[[i]] <- data.frame(levels(spectra$genotype))
  colnames(spectra.blues.list[[i]]) <- 'genotype'
  
}
#View(spectra.blues.list)

### As before, names() can be used to track the appropriate Ns.

names(spectra.blues.list) <- c('HN', 'LN')
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
    
    ### The BLUE model is 
    
  spectrum.blue.mod <- lmer(reflectance ~ genotype + (1|ASD) + (1|Rep) ,  data = temp)
    
  
     ### The BLUEs centered around the mean can also be calculated and stored in a dataframe. 
    spectra.int <- fixef(spectrum.blue.mod)[1]
    spectra.blues.temp <- fixef(spectrum.blue.mod)
    spectra.blues.temp[-1] <- spectra.blues.temp[-1] + spectra.int
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

names(spectra.blues.list[['HN']]) <- sub('.x', '', names(spectra.blues.list[['HN']]))

spectra_one_rep <- spectra[spectra$Rep == 1,]
spectra_columns <- subset(spectra_one_rep, select = c(1:12))
merged_1 <- merge(spectra_columns[which(spectra_columns$Trt== 'HN'),], spectra.blues.list[['HN']], by = 'genotype')
merged_2 <- merge(spectra_columns[which(spectra_columns$Trt== 'LN'),], spectra.blues.list[['LN']], by = 'genotype')


blues_merged <- merged_1 %>% full_join(merged_2)
blues_merged$Block <- factor(blues_merged$Block, levels= c('2', '4', '1', '3'))
write.csv(blues_merged, './spectra_blues.csv', row.names = FALSE)

blues <- read_csv('spectra_blues.csv')
blues <- as.data.frame(blues)
blues <- blues %>% select(-contains(c('.x', '.y', '.z')))

blues <- blues[ , colSums(is.na(blues))==0]


str(blues)


blues$genotype <- as.factor(blues$genotype)
blues$PLOT.ID <- as.factor(blues$PLOT.ID)
blues$rows <- as.factor(blues$rows)
blues$ranges <- as.factor(blues$ranges)
blues$Block <- as.factor(blues$Block)
blues$Rep <- as.factor(blues$Rep)
blues$Trt <- as.factor(blues$Trt)
blues$year <- as.factor(blues$year)
blues$note <-as.factor(blues$note)
blues$Calibration <- as.factor(blues$Calibration)
blues$ASD <- as.factor(blues$ASD)


blues_plot <- melt(blues, id.vars = c('PLOT.ID', 'genotype','rows', 'ranges', 'Block' ,'Rep', 'Trt', 'year', 'note', 'Calibration', 'ASD', 'Group'))
a <- sapply(strsplit(as.character(blues_plot$variable), '[.]'), '[[', 1)
blues_plot$wavelength <- as.numeric(substr(a,2,5))
blues_plot$wavelength <- as.numeric(blues_plot$wavelength)


data <- blues_plot %>% group_by(wavelength,note) %>% 
  summarise(mean.ref = mean(value, na.rm=TRUE), sd.ref = sd(value, na.rm = TRUE), se.ref= sd(value, na.rm=TRUE)/sqrt(length(value)), 
            max = max(value, na.rm = TRUE), min = min(value, na.rm = TRUE))

plt_blues <- ggplot(data=data, aes(x=wavelength, group= note)) +
  geom_line(aes(y=mean.ref, color=note), size = 0.6)+
  geom_ribbon(aes(ymin=mean.ref-se.ref , ymax=mean.ref+se.ref , fill=note),alpha=0.3)+
  labs(title = 'Leaf Spectra inbred vs hybrids under HN', caption = '**Envelopes represent 1 sd from the mean')+
  theme_bw()


plt_blues



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



