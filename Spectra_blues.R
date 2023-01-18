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
    
    spectra.blues.list[[i]] <- merge(spectra.blues.list[[i]], spectra.blues.temp, by = 'genotype', all.x= TRUE)
  
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


blues_plot <- melt(blues, id.vars = c('PLOT.ID', 'genotype','rows', 'ranges', 'Block' ,'Rep', 'Trt', 'year', 'note', 'Group'))
a <- sapply(strsplit(as.character(blues_plot$variable), '[.]'), '[[', 1)
blues_plot$wavelength <- as.numeric(substr(a,2,5))
blues_plot$wavelength <- as.numeric(blues_plot$wavelength)


data <- blues_plot %>% group_by(wavelength,Group) %>% 
  summarise(mean.ref = mean(value, na.rm=TRUE), sd.ref = sd(value, na.rm = TRUE), se.ref= sd(value, na.rm=TRUE)/sqrt(length(value)), 
            max = max(value, na.rm = TRUE), min = min(value, na.rm = TRUE))

plt_blues <- ggplot(data=data, aes(x=wavelength, group= genotype)) +
  geom_line(aes(y=mean.ref, color=Group), size = 0.6)+
  geom_ribbon(aes(ymin=mean.ref-se.ref , ymax=mean.ref+se.ref , fill=Group),alpha=0.3)+
  labs(title = 'Leaf Spectra inbred vs hybrids under HN', caption = '**Envelopes represent 1 sd from the mean')+
  theme_bw()


plt_blues

blues_plot_2 <- ggplot() +
  geom_line(data=blues_plot[which(blues_plot$Trt == 'HN'),],aes(wavelength, value, group= genotype, color='HN'), size=0.05)+
  geom_line(data=blues_plot[which(blues_plot$Trt == 'LN'),], aes(wavelength, value, group= genotype, color='LN'), size=0.05)+
  scale_color_manual(name='Treatment',  values = c('HN'= 'green', 'LN' = 'yellow'))

blues_plot_2

blues_plot <- ggplot(data = blues_plot)+
  geom_line(aes(wavelength, value, group=genotype, color=Trt), size=0.1)+
  labs(title = 'Leaf spectra (BLUES)', x='bands', y='reflectance')

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




                                  ######################  Blues for CHL #######################33


library(lme4)
library(reshape2)
library(ggplot2)
library(rrBLUP)

### The yield data can be found in the file "PhenomeForce_Yield.csv". This can be read 
### into R using read.csv().

spectra <- read.csv("Raw_spectrum_merged_predicted_CHL.csv")
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

spectra <- spectra %>% relocate(CHL , .before = 'X350')

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

View(spectra_comb)

levels(spectra_comb$Trt)

spectra_comb.list<- spectra

View(spectra_comb.list)

data<- spectra_comb.list

data.blue.mod<-lmer(CHL~ genotype + (1|Trt) + (1|genotype:Trt) + (1|Rep) , data=data)

### The fixef() function returns the fixed effects of the model. Notice that the first 
### fixed effect is the intercept. This is also the fixed effect for the first "gid" level,
### which is 4755014. 

fixef(data.blue.mod)[1:10]

data.int<-fixef(data.blue.mod)[1]
data.int

### The intercept can be added to the BLUE values to center the BLUEs around the mean

data.blues<-fixef(data.blue.mod)
data.blues[-1]<-data.blues[-1]+data.int

############# BURASI YENI EKLENDI ONEMLI !!!!!!! #########

names(data.blues)[1]<-levels(spectra$genotype)[1]
###########################################################
### The gsub() function can be used to take off the "new_GID" character string in the names
### of the BLUEs. 

names(data.blues)[1:10]
names(data.blues)[1]<-levels(data$genoype)[1]
gsub("genotype", "", names(data.blues))[1:10]

names(data.blues)<-gsub("genotype", "", names(data.blues))
data.blues[1:10]

### Using dataframe(), the data can be stored in a data frame with two columns
### for "genotype" and "grain.data.blue".

data.blues<-data.frame(names(data.blues), data.blues)
colnames(data.blues)<-c("genotype", "CHL")
head(data.blues)
rownames(data.blues) <- NULL

hist(data.blues$CHL)

write.csv(data.blues, './blues_CHL.csv', row.names = FALSE)


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
              ######################  Blues for Extra Phenos for seperate N treatments #######################



library(lme4)
library(reshape2)
library(ggplot2)
library(rrBLUP)

spectra <- read.csv("BGEM_extra_phenotypes.csv")

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

spectra <- subset(spectra, select = -c( X))

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

View(spectra)

## removing block 4  ##
spectra[c(1051:1400), c('leaf_length', 'leaf_width' , 'ear_height', 'flag_leaf',  'plant_height') ] <- NA

bands <- colnames(spectra)[c(13:17)]
bands

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
### Another list with a length of 2 can be created to store the BLUEs for the bands.

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

for(i in 2:2){  ####!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! HN IS REMOVED BECAUSE THE BLOC 4 IS MISSING !!!!!!!!!!!!!!!!!!!!!!! 
  for(j in 1:length(bands)){
    ### The hyperspectral data from just one N can be stored in a temporary variable.
    
    temp <- spectra.list[[i]]
    
    ### The hyperspectral data from just one band can be pulled out along with the field
    ### design variables. 
    
    temp <- temp[, which(colnames(temp) %in% c('genotype', 'Rep', 'Block', 'ASD', bands[j]))]
    colnames(temp)[5] <- 'reflectance'
    
    ### The BLUE model is 
    
    spectrum.blue.mod <- lmer(reflectance ~ genotype  +  (1|Rep) ,  data = temp)
    
    
    ### The BLUEs centered around the mean can also be calculated and stored in a dataframe. 
    spectra.int <- fixef(spectrum.blue.mod)[1]
    spectra.blues.temp <- fixef(spectrum.blue.mod)
    spectra.blues.temp[-1] <- spectra.blues.temp[-1] + spectra.int
    
    ############# BURASI YENI EKLENDI ONEMLI !!!!!!! #########
    names(spectra.blues.temp)[1]<-levels(spectra$genotype)[1]
    ###########################################################
    spectra.blues.temp <- data.frame (spectra.blues.temp)
    spectra.blues.temp<- cbind(genotype=rownames(spectra.blues.temp), spectra.blues.temp)
    rownames(spectra.blues.temp) <- NULL
    spectra.blues.temp$genotype <- gsub('genotype', '', spectra.blues.temp$genotype)
    colnames(spectra.blues.temp) <- c('genotype', bands[j])
    
    ### The BLUPs can be merged based on the "genotype" column with the existing dataframes 
    ### in the BLUPs list using the merge() function.
    
    spectra.blues.list[[i]] <- merge(spectra.blues.list[[i]], spectra.blues.temp, by = 'genotype', all.x = TRUE)
    
    ### counters can be used to track the progress
    
    print(i)
    print(j)
    
  }
}


names(spectra.blues.list[['HN']]) <- sub('.x', '', names(spectra.blues.list[['HN']]))

spectra_one_rep <- spectra[spectra$Rep == 1,]
spectra_columns <- subset(spectra_one_rep, select = c(1:12))
#merged_1 <- merge(spectra_columns[which(spectra_columns$Trt== 'HN'),], spectra.blues.list[['HN']], by = 'genotype', all.x = TRUE )

merged_1 <- spectra_one_rep[spectra_one_rep$Trt == 'HN', ]
merged_1 <- subset(merged_1 , select = -c(Row, Ran., Leaf.Length.1, Leaf.Width.1, Leaf.Length.2, Leaf.Width.2, Ear.Height.1, Flag.Leaf.Height.1, Plant.Height.1,
                                          Ear.Height.2, Flag.Leaf.Height.2, Plant.Height.2, Yang.Field.Notes, Inbred.or.Hybrid., Data.missing.))
merged_2 <- merge(spectra_columns[which(spectra_columns$Trt== 'LN'),], spectra.blues.list[['LN']], by = 'genotype', all.x = TRUE )


blues_merged <- merged_1 %>% full_join(merged_2)
blues_merged$Block <- factor(blues_merged$Block, levels= c('2', '4', '1', '3'))
write.csv(blues_merged, './extra_pheno_blues.csv', row.names = FALSE)

blues <- read.csv("extra_pheno_blues.csv")
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


blues_plot <- melt(blues, id.vars = c('PLOT.ID', 'genotype','rows', 'ranges', 'Block' ,'Rep', 'Trt', 'year', 'note', 'Group'))
a <- sapply(strsplit(as.character(blues_plot$variable), '[.]'), '[[', 1)
blues_plot$wavelength <- as.numeric(substr(a,2,5))
blues_plot$wavelength <- as.numeric(blues_plot$wavelength)


data <- blues_plot %>% group_by(wavelength,Group) %>% 
  summarise(mean.ref = mean(value, na.rm=TRUE), sd.ref = sd(value, na.rm = TRUE), se.ref= sd(value, na.rm=TRUE)/sqrt(length(value)), 
            max = max(value, na.rm = TRUE), min = min(value, na.rm = TRUE))

plt_blues <- ggplot(data=data, aes(x=wavelength, group= genotype)) +
  geom_line(aes(y=mean.ref, color=Group), size = 0.6)+
  geom_ribbon(aes(ymin=mean.ref-se.ref , ymax=mean.ref+se.ref , fill=Group),alpha=0.3)+
  labs(title = 'Leaf Spectra inbred vs hybrids under HN', caption = '**Envelopes represent 1 sd from the mean')+
  theme_bw()


plt_blues

blues_plot_2 <- ggplot() +
  geom_line(data=blues_plot[which(blues_plot$Trt == 'HN'),],aes(wavelength, value, group= genotype, color='HN'), size=0.05)+
  geom_line(data=blues_plot[which(blues_plot$Trt == 'LN'),], aes(wavelength, value, group= genotype, color='LN'), size=0.05)+
  scale_color_manual(name='Treatment',  values = c('HN'= 'green', 'LN' = 'yellow'))

blues_plot_2

blues_plot <- ggplot(data = blues_plot)+
  geom_line(aes(wavelength, value, group=genotype, color=Trt), size=0.1)+
  labs(title = 'Leaf spectra (BLUES)', x='bands', y='reflectance')

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





######################  Blues for CHL for seperate N TReatmens #######################33


library(lme4)
library(reshape2)
library(ggplot2)
library(rrBLUP)

### The yield data can be found in the file "PhenomeForce_Yield.csv". This can be read 
### into R using read.csv().

spectra <- read.csv("Raw_CHL_predicted")


str(spectra)

spectra$Rep <- factor(spectra$Rep)
spectra$genotype <- factor(spectra$genotype)
spectra$Trt <- factor(spectra$Trt)
spectra$Group <- factor(spectra$Group)
spectra$PLOT.ID <- factor(spectra$PLOT.ID)


View(spectra)



spectra_comb <- spectra[spectra$Trt == 'LN' & spectra$Group == 'Inbred' ,]

levels(spectra_comb$Trt)

spectra_comb.list<- spectra_comb

data<- spectra_comb.list

data.blue.mod<-lmer(CHL~ genotype + (1|Rep) + (1|Rep:genotype) , data=data)

### The fixef() function returns the fixed effects of the model. Notice that the first 
### fixed effect is the intercept. This is also the fixed effect for the first "gid" level,
### which is 4755014. 

fixef(data.blue.mod)[1:10]

data.int<-fixef(data.blue.mod)[1]
data.int

### The intercept can be added to the BLUE values to center the BLUEs around the mean

data.blues<-fixef(data.blue.mod)
data.blues[-1]<-data.blues[-1]+data.int

############# BURASI YENI EKLENDI ONEMLI !!!!!!! #########

names(data.blues)[1]<-levels(spectra$genotype)[1]
###########################################################
### The gsub() function can be used to take off the "new_GID" character string in the names
### of the BLUEs. 

names(data.blues)[1:10]
gsub("genotype", "", names(data.blues))[1:10]

names(data.blues)<-gsub("genotype", "", names(data.blues))
data.blues[1:10]

### Using dataframe(), the data can be stored in a data frame with two columns
### for "genotype" and "grain.data.blue".

data.blues<-data.frame(names(data.blues), data.blues)
colnames(data.blues)<-c("genotype", "CHL")
head(data.blues)
rownames(data.blues) <- NULL

hist(data.blues$CHL)

write.csv(data.blues, './blues_CHL_HN_only_inbreds.csv', row.names = FALSE)
write.csv(data.blues, './blues_CHL_LN_only_inbreds.csv', row.names = FALSE)

ggplot(data=spectra)+
  geom_density(aes(CHL, fill=Group), alpha=0.3)+
  theme_classic()+
  labs(title='CHL contents Inbreds vs Hybrids')+
  scale_fill_brewer(palette = 'Spectral')


data_LN <- read.csv('blues_CHL_LN_only_inbreds.csv')
data_HN <- read.csv('blues_CHL_HN_only_inbreds.csv') 

names(data_LN)[2] <- 'CHL_LN'
names(data_HN)[2] <- 'CHL_HN'


blues_CHL_LN_HN_sperately <- merge(data_LN, data_HN , by = 'genotype')

View(blues_CHL_LN_HN_sperately)

write.csv(blues_CHL_LN_HN_sperately, './blues_CHL_LN_HN_seperately.csv')
