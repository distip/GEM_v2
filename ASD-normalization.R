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
