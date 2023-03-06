library(Matrix)
library(lme4)
library(reshape2)
library(rrBLUP)
library(ggplot2)
library(tidyverse)
library(viridisLite)
library(viridis)
library(grid)
library(lmerTest)

setwd('/home/schnablelab/deniz/GEM_v2/GEM_v2')
data <- read_csv('Spectrum_with_biochemical_all_predicted.csv')

View(data)

colnames(data)

data = data[, c("PLOT ID" , "Group" , "genotype",  "Trt",         "% N",         "% P",         "% K" ,        "% S",         "% Ca",       
                "% Mg" ,       "ppm Zn" ,     "ppm Fe" ,     "ppm Mn"  ,    "ppm Cu"  ,    "ppm B" ,      "ppm Mo" )]

data_plot <- data %>% melt(id.vars= c('PLOT ID', 'Group', 'genotype', 'Trt'))
View(data_plot)
colnames(data_plot)[5:6] <- c('Trait', 'value')

ggplot(data_plot[data_plot$Group == 'Hybrid', ])+
  geom_density(aes(x=value, fill=Trt), alpha=0.5)+
  facet_wrap('Trait', scales = 'free')

ggplot(data_plot[data_plot$Group == 'Inbred', ])+
  geom_boxplot(aes(y=value, fill=Trt), alpha=0.5)+
  facet_wrap('Trait', scales = 'free')



                                  ############## BLUPS ###############

spectra <- read.csv("Spectrum_with_biochemical_all_predicted.csv")
spectra <- spectra#[spectra$Group == 'Hybrid', ]
spectra <- spectra %>% select(genotype, Trt, Group, Rep, X..N, X..P, X..K, X..S, X..Ca, 
                              X..Mg, ppm.Zn, ppm.Fe, ppm.Mn, ppm.Cu, ppm.B, ppm.Mo)

str(spectra)

spectra$Rep <- factor(spectra$Rep)
spectra$genotype <- factor(spectra$genotype)
spectra$Trt <- factor(spectra$Trt)
spectra$Group <- factor(spectra$Group)
#spectra$PLOT.ID <- factor(spectra$PLOT.ID)

#View(spectra)

#levels(spectra$Trt)

spectra.list <- vector('list' , 2)
#spectra.list

for(i in 1:2) {
  ##The appropriate N can be pulled out using which().
  
  spectra.list[[i]] <- spectra[which(spectra$Trt ==levels(spectra$Trt)[i]),]
  
  ### To keep track of where each N appears within the list, names() can be used to 
  ### set the name of each list element to the N
  
  names(spectra.list)[i] <- levels(spectra$Trt)[i]
  
}

#View(spectra.list)

### A dataframe can be created to store the broad-sense heritability estimates for each
### band on each N. The expand.grid() function creates a data frame where each row 
### contains a different combination of the bands and N. An "H2" column can be added
### to store the broad-sense heritability estimates. 

bands <- colnames(spectra)[5:16]
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
#View(spectra.blups.list)

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
    
    temp <- temp[, which(colnames(temp) %in% c('genotype', 'Rep', bands[j]))]
    colnames(temp)[3] <- 'reflectance'
    
    ### The BLUP model is 
    
    spectrum.blup.mod <- lmer(reflectance ~ (1|genotype) +(1|Rep) , data = temp)
    
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
View(spectra.blups.list)
head(bands.H2)

### To plot different lines for each of the different N application, the color can be set to "Trt".

ggplot(bands.H2, aes(y=H2, x=Trt)) + geom_bar(stat='identity')+
  labs(title = 'Heritability of  (Only Hybrids)')+
  ylim(c(0, 0.8))+
  theme_bw(10)+
  facet_wrap('bands', scales='free')

CHL_blups_HN <- spectra.blups.list[['HN']]
CHL_blups_LN <- spectra.blups.list[['LN']]
names(CHL_blups_HN)[2] <- 'CHL_blups_HN'
names(CHL_blups_LN)[2] <- 'CHL_blups_LN'

blups_CHL_LN_HN_sperately <- merge(CHL_blups_HN, CHL_blups_LN , by = 'genotype')
View(blups_CHL_LN_HN_sperately)

write.csv(blups_CHL_LN_HN_sperately, './blups_CHL_LN_HN_seperately.csv', row.names = FALSE)


aa <- blups_CHL_LN_HN_sperately %>% melt(id.vars = 'genotype')

aa %>% ggplot(aes(y=value, x=variable))+
  geom_boxplot(aes(fill=variable))

