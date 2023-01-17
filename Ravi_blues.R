

setwd('/home/schnable/Documents/Ravi_blues/')

library(lme4)
library(reshape2)
library(ggplot2)
library(rrBLUP)

### The yield data can be found in the file "PhenomeForce_Yield.csv". This can be read 
### into R using read.csv().

spectra <- read.csv("Bugeater2021_merged_v6.2.csv")
View(spectra)

str(spectra)

spectra$Rep <- factor(spectra$Rep)
spectra$Row <- factor(spectra$Row)
spectra$Genotype <- factor(spectra$Genotype)
spectra$Column <- factor(spectra$Column)



bands <- colnames(spectra)[c(7:41)]
bands <- bands[-3]

spectra.list <- vector('list' , 1)
spectra.list

for(i in 1:1) {
  ##The appropriate N can be pulled out using which().
  
  spectra.list[[i]] <- spectra
  
  ### To keep track of where each N appears within the list, names() can be used to 
  ### set the name of each list element to the N
  
  names(spectra.list)[i] <- 'blues'
  
}


View(spectra.list)

spectra.blues.list <- vector('list', 1)
spectra.blues.list
names(spectra.blues.list) <- 'blues'

for(i in 1:1){
  
  ### To get started, a dataframe containing only the "gid" column can be created for each
  ### element of the list. 
  
  spectra.blues.list[[i]] <- data.frame(levels(spectra$Genotype))
  colnames(spectra.blues.list[[i]]) <- 'Genotype'
  
}

view(spectra.blues.list)
### The following nested loop will first go through the list of dataframes each containing
### hyperspectral data from a single N.

for(i in 1:1){
  for(j in 1:length(bands)){
    ### The hyperspectral data from just one N can be stored in a temporary variable.
    
    temp <- spectra.list[[i]]
    
    ### The hyperspectral data from just one band can be pulled out along with the field
    ### design variables. 
    
    temp <- temp[, which(colnames(temp) %in% c('Genotype', 'Rep',  bands[j]))]
    colnames(temp)[3] <- 'reflectance'
    
    ### The BLUE model is 
    
    spectrum.blue.mod <- lmer(reflectance ~ Genotype  +  (1|Rep) ,  data = temp)
    
    
    ### The BLUEs centered around the mean can also be calculated and stored in a dataframe. 
    spectra.int <- fixef(spectrum.blue.mod)[1]
    spectra.blues.temp <- fixef(spectrum.blue.mod)
    spectra.blues.temp[-1] <- spectra.blues.temp[-1] + spectra.int
    
    ############# BURASI YENI EKLENDI ONEMLI !!!!!!! #########
    names(spectra.blues.temp)[1]<-levels(spectra$Genotype)[1]
    ###########################################################
    spectra.blues.temp <- data.frame (spectra.blues.temp)
    spectra.blues.temp<- cbind(Genotype=rownames(spectra.blues.temp), spectra.blues.temp)
    rownames(spectra.blues.temp) <- NULL
    spectra.blues.temp$genotype <- gsub('Genotype', '', spectra.blues.temp$Genotype)
    colnames(spectra.blues.temp) <- c('Genotype', bands[j])
    
    ### The BLUPs can be merged based on the "genotype" column with the existing dataframes 
    ### in the BLUPs list using the merge() function.
    
    spectra.blues.list[[i]] <- merge(spectra.blues.list[[i]], spectra.blues.temp, by = 'Genotype')
    
    ### counters can be used to track the progress
    
    print(i)
    print(j)
    
  }
}

view(spectra.blues.list)
