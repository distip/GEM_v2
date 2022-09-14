library(Matrix)
library(lme4)
library(reshape2)
library(rrBLUP)
library(ggplot2)
library(glmmTMB)
library(ncf)
library(viridisLite)
library(viridis)
library(tibble)
library(PerformanceAnalytics)
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
spectra$rows <- as.numeric(spectra$rows)
spectra$ranges <- as.numeric(spectra$ranges)
spectra$PLOT.ID <- factor(spectra$PLOT.ID)

spectra <- subset(spectra, select = -c(X, Unnamed..0))


spectra.new <- spectra
spectra.new$Block <- factor(spectra.new$Block, levels= c('2', '4', '1', '3'))

                                                 ### heat map of the field for one wavelngth ###

ggplot(spectra.new, aes(rows, ranges, color=X730)) + 
  geom_point(size=1.3) +
  scale_y_continuous(name='ranges', limits = c(1,13))+
  geom_rect( aes(xmin=0.4, xmax = 51, ymin = 7.5, ymax = 12.5), fill=NA, colour='red')+
  #annotate('rect', xmin=0, xmax = 50, ymin = 7.5, ymax = 12.5, alpha= .1)+
  facet_wrap(Block~Trt)+
  scale_colour_viridis() +
  labs(title = 'Range and Row effects on leaf spectrum', caption = ' Blocks 1 and 4 = + N , 2 and 3 = -N\nred rectangles are the hybrids ')+
  theme_classic()
  #theme(strip.background = element_rect(color = 'black', fill = Trt))


                                ######## Calculating the spatial effects (do this before blups calculation) ###########

spectra <- add_column(spectra, pos = numFactor(scale(as.numeric(spectra$rows)), scale(as.numeric(spectra$ranges))), .after = 'ranges' )
spectra <- add_column(spectra, ID = factor(rep(1, nrow(spectra))), .after= 'pos')

data <- spectra[which(spectra$Block == 4),]
data <- spectra
m1 <-glmmTMB(X730 ~ (1|genotype), data= data)
m2 <-glmmTMB(X730 ~ (1|genotype) + ar1(pos + 0 | ID) , data= data)
m3 <-glmmTMB(X730 ~ (1|genotype) + (1|ASD) + ar1(pos + 0 | ID), data=data)


anova(m1,m2,m3)

                            ####################### Calculation of BLUPs form the m2 (only spatial) for each blocks ########################

levels(spectra$Block)

spectra.list <- vector('list' , 4)
spectra.list

for(i in 1:4) {
  ##The appropriate N can be pulled out using which().
  
  spectra.list[[i]] <- spectra[which(spectra$Block ==levels(spectra$Block)[i]),]
  
  ### To keep track of where each N appears within the list, names() can be used to 
  ### set the name of each list element to the N
  
  names(spectra.list)[i] <- levels(spectra$Block)[i]
  
}

View(spectra.list)

### A dataframe can be created to store the broad-sense heritability estimates for each
### band on each N. The expand.grid() function creates a data frame where each row 
### contains a different combination of the bands and N. An "H2" column can be added
### to store the broad-sense heritability estimates. 

bands <- colnames(spectra)[-c(1:13)]
bands

Block <- as.character(levels(spectra$Block))
Block

bands.H2 <- expand.grid(bands, Block)
bands.H2

colnames(bands.H2) <- c('bands', 'Block')
bands.H2$H2 <- NA
head(bands.H2)

### Another list with a length of 2 can be created to store the BLUPs for the bands.

spectra.blups.list <- vector('list', 4)
spectra.blups.list

for(i in 1:4){
  
  ### To get started, a dataframe containing only the "genotype" column can be created for each
  ### element of the list. 
  
  spectra.blups.list[[i]] <- data.frame(levels(spectra$genotype))
  colnames(spectra.blups.list[[i]]) <- 'genotype'
  
}
View(spectra.blups.list)

### As before, names() can be used to track the appropriate Ns.

names(spectra.blups.list) <- c('1', '2', '3', '4')
### The following nested loop will first go through the list of dataframes each containing
### hyperspectral data from a single N.

for(i in 1:length(spectra.list)){
  for(j in 1:length(bands)){
    ### The hyperspectral data from just one block can be stored in a temporary variable.
    
    temp <- spectra.list[[i]]
    
    ### The hyperspectral data from just one band can be pulled out along with the field
    ### design variables. 
    
    temp <- temp[, which(colnames(temp) %in% c('genotype', 'pos', 'ID', 'ASD', 'Block', bands[j]))]
    colnames(temp)[6] <- 'reflectance'
    
    ### The BLUP model is 
    
    spectrum.blup.mod <- glmmTMB(reflectance ~ (1|genotype) + ar1(pos + 0 | ID) , data= temp)
    
    
    ### The variance components can be extracted to calculate broad-sense heritability.
    
    #Vg <- data.frame(VarCorr(spectrum.blup.mod))$vcov[1]
    #Ve <- data.frame(VarCorr(spectrum.blup.mod))$vcov[3]
    #H2 <- Vg/(Vg+(Ve/2))
    
    
    ### The broad-sense heritability can be put in the appropriate place within the 
    ### bands.H2 dataframe using which().
    
    #bands.H2[which(bands.H2$bands == bands[j] & bands.H2$Block == names(spectra.list)[i]), 'H2'] <- H2
    
    ### The BLUPs centered around the mean can also be calculated and stored in a dataframe. 
    
    spectra.blups.temp <- ranef(spectrum.blup.mod)$cond$genotype + fixef(spectrum.blup.mod)$cond
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




### To look at the hyperspectral data from multiple genotypes on a single data, the first 
### fifteen genotypes from the first date can be pulled out. As before, melt() moves the 
### spectral bands from columns to rows, and the wavelengths need to be extracted from the
### band labels.


spectra.blups.sub.melt <- melt(spectra.blups.list)
#a <- sapply(strsplit(as.character(spectra.blups.sub.melt$variable), '[.]'), '[[', 1)
spectra.blups.sub.melt$wavelength <- as.numeric(substr(spectra.blups.sub.melt$variable,2,5))



View(spectra.blups.sub.melt)
View(spectra.blups.list)
colnames(spectra.blups.sub.melt)[4] <- 'Block'


write.csv(spectra.blups.sub.melt, './fullspectra_plots_spatial_cor.csv', row.names = FALSE )
spectra.blups.sub.melt <- read.csv('fullspectra_plots.csv')


spectra_columns <- subset(spectra, select = c(1:13))
merged_1 <- merge(spectra_columns[which(spectra_columns$Block== '1'),], spectra.blups.list[['1']], by = 'genotype')
merged_2 <- merge(spectra_columns[which(spectra_columns$Block== '2'),], spectra.blups.list[['2']], by = 'genotype')
merged_3 <- merge(spectra_columns[which(spectra_columns$Block== '3'),], spectra.blups.list[['3']], by = 'genotype')
merged_4 <- merge(spectra_columns[which(spectra_columns$Block== '4'),], spectra.blups.list[['4']], by = 'genotype')

blups_merged <- merged_1 %>% full_join(merged_2) %>% full_join(merged_3) %>% full_join(merged_4) 
blups_merged$Block <- factor(blups_merged$Block, levels= c('2', '4', '1', '3'))

ggplot(blups_merged, aes(rows, ranges, color=X730)) + 
  geom_point(size=1.3) +
  scale_y_continuous(name='ranges', limits = c(1,13))+
  geom_rect( aes(xmin=0.4, xmax = 51, ymin = 7.5, ymax = 12.5), fill=NA, colour='red')+
  #annotate('rect', xmin=0, xmax = 50, ymin = 7.5, ymax = 12.5, alpha= .1)+
  facet_wrap(Block~Trt)+
  scale_colour_viridis() +
  labs(title = 'Range and Row effects on leaf spectrum after spatial correction', caption = ' Blocks 1 and 4 = + N , 2 and 3 = -N\nred rectangles are the hybrids ')+
  theme_classic()
#theme(strip.background = element_rect(color = 'black', fill = Trt))




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




                                  ####################### Variance Partitioning ########################

library(Matrix)
library(lme4)
library(reshape2)
library(rrBLUP)
library(ggplot2)
library(tidyverse)



spectra <- blups_merged
str(spectra)

bands <- (colnames(spectra)[-c(1:13)])
bands_df <-as.data.frame(colnames(spectra)[-c(1:13)])
View(bands_df) 
bands_df$H2 <- NA
colnames(bands_df) <- c('band', 'H2')
bands.H2 <- bands_df

spectral.blups.list <- data.frame(levels(spectra$genotype))

colnames(spectral.blups.list) <- 'genotype'


extractVarsLmm <- function(model, id_genotype = c("genotype"), 
                           overdisp = "PLOT.ID"){
  
  # Get variance covariance matrix of random terms
  vcov_all <- VarCorr(model)
  
  # Get the variance-covariance matrix from line ID
  vcov_matrix <- vcov_all[[which(names(vcov_all) %in% id_genotype)]] %>% 
    as.data.frame %>% as.matrix()
  
  
  #### Calculate all variance-covariance components ####
  # Get the variance in slopes (plasticity) from the model output
  var_slopes <- vcov_matrix[2, 2]
  
  # Calculate the variance on each nitrate and covariance between them
  ## equation 4 in https://link.springer.com/article/10.1007/s00265-013-1603-9
  phi <- matrix(c(1, 1, 0, 1), nrow = 2) # design matrix
  env_vcov_matrix <- phi %*% vcov_matrix %*% t(phi)
  
  var_ln <- env_vcov_matrix[1, 1] 
  var_hn <- env_vcov_matrix[2, 2]
  cov_ln_hn <- env_vcov_matrix[1, 2]
  
  
  # Calculate covariance between each nitrate and the slopes 
  ## note that: 
  ## COV(LN; slopes) = COV(HN; LN) - V(LN)
  ## and equally:
  ## COV(HN; slopes) = COV(HN; LN) - V(HN)
  ## This follows from looking at the matrix calculation in equation 4 of Bromer 2013. 
  ## Also, it matches output from lmer, for example looking at vcov_matrix[1,2], which gives COV(LN; slopes)
  cov_ln_slopes <- cov_ln_hn - var_ln
  cov_hn_slopes <- cov_ln_hn - var_hn
  
  
  #### Calculate correlations ####
  # From all these variances and covariances we can calculate several correlations
  ## remembering that COR12 = COV12/sqrt(V1 * V2)
  
  # across-environment correlation
  cor_ln_hn <- cov_ln_hn/sqrt(var_ln * var_hn)
  
  # correlation between environemnt and slopes
  cor_ln_slopes <- (cov_ln_slopes)/sqrt(var_ln * var_slopes)
  cor_hn_slopes <- (cov_hn_slopes)/sqrt(var_hn * var_slopes)
  
  #### Calculate residual variance ####
  # Extract the residual variance
  if(family(model)$family == "gaussian"){
    var_res <- attr(VarCorr(model), "sc")^2
  } else if(family(model)$family == "poisson" & family(model)$link == "log"){
    
    # Calculate distribution-specific variance following table 2 of Nakagawa & Schielzeth 2010
    # Beta0 adapted from the function `r.squaredGLMM.merMod` from the MuMIn package
    beta0 <- mean(model.matrix(model) %*% fixef(model))
    var_dis <- log(1/exp(beta0)+1)
    
    # Add the over-dispersion variance to get total "residual" variance
    var_res <- var_dis + vcov_all[[overdisp]][1]
  }
  
  #### Calculate fixed effects variance ####
  var_fixed <- var(model.matrix(model) %*% fixef(model))
  
  
  
  #### Tidy output #####
  # Output a data.frame
  data.frame(var_ln = var_ln,
             var_hn = var_hn,
             var_nitrate = var_fixed,
             var_plas = var_slopes,
             var_res = var_res)
}



var.part.list <- vector('list', length(bands))
names(var.part.list) <- as.numeric(substr(as.character(bands), 2,5))

values <- c()
for(i in 1:length(bands)){
  temp <- spectra
  temp<-temp[, which(colnames(spectra) %in% c('PLOT.ID','genotype', 'Block', 'Trt', 'Rep', 'ASD', bands[i]))]
  colnames(temp)[7] <- 'reflectance'
  m <- lmer(reflectance ~ Trt + (Trt|genotype), data=temp)
  extractVarsLmm(m)
  LN <- extractVarsLmm(m)[[1]]/ sum(extractVarsLmm(m))
  HN <- extractVarsLmm(m)[[2]]/ sum(extractVarsLmm(m))
  Nitrate <- extractVarsLmm(m)[[3]]/ sum(extractVarsLmm(m))
  Plasticity <- extractVarsLmm(m)[[4]]/ sum(extractVarsLmm(m))
  Res <- extractVarsLmm(m)[[5]]/ sum(extractVarsLmm(m))
  values <-  c(LN, HN, Nitrate, Plasticity, Res)
  var.part.list[[i]] <- values
  
}


var.part.list.melt <- melt(var.part.list)
var.part.list.melt$source <- rep(c('LN', 'HN', 'Nitrate', 'Plasticity', 'Residual'), 10755/5)
var.part.list.melt$source <- factor(var.part.list.melt$source)
var.part.list.melt$band <- as.numeric(var.part.list.melt$band)
colnames(var.part.list.melt) <- c('values', 'band', 'source')


ggplot(var.part.list.melt, aes(fill= source, y=values, x=band)) +
  geom_bar(position='fill', stat='identity', width = 1, alpha=0.8)+
  theme_classic()+
  theme(legend.background = element_rect(fill = '#FFCC66',color='grey50',  size=1))+
  scale_fill_brewer(palette= 'Set3')+
  labs(title = 'Variance Partitioning of Leaf Spectrum', y ='%', x='wavelengths')



                                            #######Calculation of heritabilities##############


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

bands <- colnames(spectra)[-c(1:13)]
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
    
    spectrum.blup.mod <- lmer(reflectance ~ (1|genotype)+(1|ASD)+(1|Block)+(1|Rep:Block), data = temp)
    
    ### The variance components can be extracted to calculate broad-sense heritability.
    
    Vg <- data.frame(VarCorr(spectrum.blup.mod))$vcov[1]
    Ve <- data.frame(VarCorr(spectrum.blup.mod))$vcov[5]
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
  labs(title = 'Broad Sense Heritabilities Under Different Nitrogen Applications after spatial correction')




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



