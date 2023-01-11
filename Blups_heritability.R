
## This will be your trait data
spectra <- read.csv("Raw_spectrum_merged.csv")

#levels of the factor that you want to monitor (in my case I calculated blups for 2 N treatments (Trt))
levels(spectra$Trt)

# Since my data set  had 2 Nitrogen treatments (trt), I created a list with the length 2. We will store our blups here.
spectra.list <- vector('list' , 2)
spectra.list

for(i in 1:2) {
  ##The appropriate treatment level(year or N treatment) can be pulled out from your main df by using which().
  
  spectra.list[[i]] <- spectra[which(spectra$Trt ==levels(spectra$Trt)[i]),]
  
  ### To keep track of where each N appears within the list, I used names() to 
  ### set the name of each list element to the N
  
  names(spectra.list)[i] <- levels(spectra$Trt)[i]
  
}

View(spectra.list)

### A dataframe can be created to store the broad-sense heritability estimates for each
### band on each N. The expand.grid() function creates a data frame where each row 
### contains a different combination of the bands (it will be traits in your case) and Nitrogen. An "H2" column can be added
### to store the broad-sense heritability estimates. 

bands <- colnames(spectra)[-c(1:11)] ## specifyy the columns of traits in your dataset
bands

Trt <- as.character(levels(spectra$Trt)) ##levels of the factor
Trt

bands.H2 <- expand.grid(bands, Trt)
bands.H2

colnames(bands.H2) <- c('bands', 'Trt')
bands.H2$H2 <- NA
head(bands.H2)

### Another list with a length of 2 (because I had two Nitrogen application) can be created to store the BLUPs for the bands (traits).

spectra.blups.list <- vector('list', 2)
spectra.blups.list

for(i in 1:2){
  
  ### To get started, a dataframe containing only the '"gid"'genotype id" column can be created for each
  ### element of the list. 
  
  spectra.blups.list[[i]] <- data.frame(levels(spectra$genotype)) ##create a dataframe in the list you created containing only the geno ids
  colnames(spectra.blups.list[[i]]) <- 'genotype' # give an appropriate name to the column (it should be the same that you have in your main dataframe)
  
}
View(spectra.blups.list)

### As before, names() can be used to track the appropriate Ns.

names(spectra.blups.list) <- c('HN', 'LN')
### The following nested loop will first go through the list of dataframes each containing
### hyperspectral data from a single N.

#This loop will calculate the blups and keep track them seperately.

for(i in 1:length(spectra.list)){
  for(j in 1:length(bands)){
    ### The hyperspectral data from just one N can be stored in a temporary variable.
    
    temp <- spectra.list[[i]]
    
    ### The hyperspectral data from just one band can be pulled out along with the field
    ### design variables. 
    
    temp <- temp[, which(colnames(temp) %in% c('genotype', 'Rep', 'Block', 'ASD', bands[j]))]
    colnames(temp)[5] <- 'reflectance'
    
    ### The BLUP model is 
    
    spectrum.blup.mod <- lmer(reflectance ~ (1|genotype)+(1|ASD)+(1|Rep) , data = temp) #you can add more effect and also interactions here such as (1|genotype : Rep)
    
    ### The variance components can be extracted to calculate broad-sense heritability.
    
    Vg <- data.frame(VarCorr(spectrum.blup.mod))$vcov[1]
    Ve <- data.frame(VarCorr(spectrum.blup.mod))$vcov[4] ## 4th elementh was ued here since I have three effects in my model
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

### To plot the broad-sense heritabilities for each band, the wavelengths need to be 
### extracted from the band labels as before.

bands.H2$bands <- as.character(bands.H2$bands)
bands.H2$bands<- as.numeric(substr(bands.H2$bands,2,5))

### To plot different lines for each of the different N application, the color can be set to "Trt". I used line plot since I have 2000 traits but you can usebarplot.

ggplot(bands.H2, aes(bands, H2)) + geom_line()+
  labs(title = 'Broad Sense Heritabilities Under Different Nitrogen Applications')+
  theme_bw(16)

#I merged the corrected blups for 2 nitrogen apps.
spectra_columns <- subset(spectra, select = c(1:11))
merged_1 <- merge(spectra_columns[which(spectra_columns$Trt== 'HN'),], spectra.blups.list[['HN']], by = 'genotype')
merged_2 <- merge(spectra_columns[which(spectra_columns$Trt== 'LN'),], spectra.blups.list[['LN']], by = 'genotype')
blups_merged <- merged_1 %>% full_join(merged_2)
write.csv(blups_merged, './spectra_blups.csv', row.names = FALSE)
