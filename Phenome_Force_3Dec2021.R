library(lme4)
library(reshape2)
library(ggplot2)
library(rrBLUP)

### The yield data can be found in the file "PhenomeForce_Yield.csv". This can be read 
### into R using read.csv().

yield<-read.csv("PhenomeForce_Yield.csv")
head(yield)

hist(yield$grain.yield)

### The yield data contains 39 trials. Each trial contains different genotypes except for 
### two repeated checks: gids 4755014 and 6176013

boxplot(yield$grain.yield~yield$trial)

table(yield$gid)[1:10]

### The variables "gid", "trial", "rep", and "block"  are currently characters/integers. 
### These need to be categorical variables to be able to analyze the field design. The 
### factor() function changes these variables to categorical variables, or "factors".

str(yield)

yield$gid<-factor(yield$gid)
yield$trial<-factor(yield$trial)
yield$rep<-factor(yield$rep)
yield$block<-factor(yield$block)

### BLUEs for grain yield can be calculated using the lmer() function of the "lme4"
### package. "gid" is set as a fixed effect, while "trial", "rep", and "block" are 
### set as random effects. "rep" is nested within "trial", and "block" is nested within
### "rep" and "trial".

yield.blue.mod<-lmer(grain.yield~gid+(1|trial)+(1|rep:trial)+(1|block:rep:trial), data=yield)

### The fixef() function returns the fixed effects of the model. Notice that the first 
### fixed effect is the intercept. This is also the fixed effect for the first "gid" level,
### which is 4755014. 

fixef(yield.blue.mod)[1:10]

yield.int<-fixef(yield.blue.mod)[1]
yield.int

### The intercept can be added to the BLUE values to center the BLUEs around the mean

yield.blues<-fixef(yield.blue.mod)
yield.blues[-1]<-yield.blues[-1]+yield.int

### The gsub() function can be used to take off the "gid" character string in the names
### of the BLUEs. 

names(yield.blues)[1:10]
names(yield.blues)[1]<-levels(yield$gid)[1]
gsub("gid", "", names(yield.blues))[1:10]

names(yield.blues)<-gsub("gid", "", names(yield.blues))
yield.blues[1:10]

### Using dataframe(), the data can be stored in a data frame with two columns
### for "gid" and "grain.yield.blue".

yield.blues<-data.frame(names(yield.blues), yield.blues)
colnames(yield.blues)<-c("gid", "grain.yield.blue")
head(yield.blues)

hist(yield.blues$grain.yield.blue)

### BLUPs can be calculated for grain yield using the same model as before but with
### "gid" changed to a random effect

yield.blup.mod<-lmer(grain.yield~(1|gid)+(1|trial)+(1|rep:trial)+(1|block:rep:trial), data=yield)

### The ranef() function returns the random effects. Since multiple variables were fit
### as random, the "$" can be used to return just the random effects for "gid".

ranef(yield.blup.mod)
ranef(yield.blup.mod)$gid

### In this model, the only fixed effect is the intercept. fixef() returns this. To
### center the BLUPs around the mean, the intercept can be added.

fixef(yield.blup.mod)

yield.blups<-ranef(yield.blup.mod)$gid+fixef(yield.blup.mod)
head(yield.blups)

### Using dataframe(), the data can be stored in a data frame with two columns
### for "gid" and "grain.yield.blup".

yield.blups<-data.frame(rownames(yield.blups), yield.blups)
colnames(yield.blups)<-c("gid", "grain.yield.blup")
head(yield.blups)

hist(yield.blups$grain.yield.blup)

### summary() shows various information about the model, including the variance 
### components for the random effects. Unfortunately, these cannot be extracted 
### directly from summary(). The VarCorr() function is needed for this but does not
### display the variance components ("vcov") unless it is inside the data.frame() 
### function.

summary(yield.blup.mod)
VarCorr(yield.blup.mod)
data.frame(VarCorr(yield.blup.mod))

### Pulling out the first (for "gid") and fifth (for error) variance components
### allows for the calculation of broad-sense heritability.

Vg<-data.frame(VarCorr(yield.blup.mod))$vcov[1]
Ve<-data.frame(VarCorr(yield.blup.mod))$vcov[5]

Vg/(Vg+(Ve/3))

### The hyperspectral data can be found in the file "PhenomeForce_Spectra.csv". This 
### can be read into R using read.csv(). The first four columns are the same as those
### found in the yield data set. The column "Date" indicates the date of phenotyping.
### Each hyperspectral band has a column labeled with the band number (e.g., "B5")
### followed by the wavelength it corresponds to (e.g., "B5.427"). There are 61 bands.

spectra<-read.csv("PhenomeForce_Spectra.csv")
View(spectra)

### As with the yield dataset, the "gid", "trial", "rep", and "block" are currently 
### characters/integers and must be changed to categorical variables. The "Date" variable 
### must also be changed to a categorical variable. 

str(spectra)

spectra$gid<-factor(spectra$gid)
spectra$trial<-factor(spectra$trial)
spectra$rep<-factor(spectra$rep)
spectra$block<-factor(spectra$block)
spectra$Date<-factor(spectra$Date)

### To look at the spectral data for a single genotype, the which() function can be used
### to pull out only observations for the gid 7305427. Using the melt() function, the 
### spectral band data will move from columns to rows.

spectra.sub<-spectra[which(spectra$gid==7305427),]
spectra.sub.melt<-melt(spectra.sub, id.vars = c("gid", "trial", "rep", "block", "Date"))
head(spectra.sub.melt)

colnames(spectra.sub.melt)[6:7]<-c("band", "reflectance")

### To plot the spectral data, the wavelength needs to be extracted from the band label
### by taking off the band number. 

spectra.sub.melt$band<-as.character(spectra.sub.melt$band)
spectra.sub.melt$band<-as.numeric(sapply(strsplit(spectra.sub.melt$band, "[.]"), "[", 2))
spectra.sub.melt$band

### Since the field design has not yet been analyzed, there are still three observations 
### for this genotype on each data (for the three replicates). To plot just the first 
### replicate, the which() function can be used to call only those data from the first 
### replicate. In aes(), the x- and y-axes can be set to "band" and "reflectance". To 
### plot the different dates as different lines, the color can be set to "Date". 

ggplot(data=spectra.sub.melt[which(spectra.sub.melt$rep==1),], aes(band, reflectance, color=Date)) + geom_line()

### To analyze the field design for the spectral dataset, it is easiest to separate the 
### different dates to create an individual dataframe for each date. These can be stored
### in a list. Since there are five dates, the list length is five. 

levels(spectra$Date)

spectra.list<-vector("list", 5)
spectra.list

for(i in 1:5){
  
  ### The appropriate date can be pulled out using which().
  
  spectra.list[[i]]<-spectra[which(spectra$Date==levels(spectra$Date)[i]),]
  
  ### To keep track of where each date appears within the list, names() can be used to 
  ### set the name of each list element to the date.
  
  names(spectra.list)[i]<-levels(spectra$Date)[i]
  
}

View(spectra.list)

### A dataframe can be created to store the broad-sense heritability estimates for each
### band on each date. The expand.grid() function creates a data frame where each row 
### contains a different combination of the bands and dates. An "H2" column can be added
### to store the broad-sense heritability estimates. 

bands<-colnames(spectra)[-c(1:5)]
bands

dates<-as.character(levels(spectra$Date))
dates

bands.H2<-expand.grid(bands, dates)
head(bands.H2)

colnames(bands.H2)<-c("band", "date")
bands.H2$H2<-NA
head(bands.H2)

### Another list with a length of five can be created to store the BLUPs for the bands. 

spectra.blups.list<-vector("list", 5)
spectra.blups.list

for(i in 1:5){
  
  ### To get started, a dataframe containing only the "gid" column can be created for each
  ### element of the list. 
  
  spectra.blups.list[[i]]<-data.frame(levels(spectra$gid))
  colnames(spectra.blups.list[[i]])<-"gid"

}

head(spectra.blups.list[[1]])

### As before, names() can be used to track the appropriate dates. 

names(spectra.blups.list)<-names(spectra.list)

### The following nested loop will first go through the list of dataframes each containing
### hyperspectral data from a single date. 

for(i in 1:length(spectra.list)){
  
  for(j in 1:length(bands)){
    
    ### The hyperspectral data from just one date can be stored in a temporary variable.
    
    temp<-spectra.list[[i]]

    ### The hyperspectral data from just one band can be pulled out along with the field
    ### design variables. 
    
    temp<-temp[,which(colnames(temp)%in%c("gid", "trial", "rep", "block", bands[j]))]
    colnames(temp)[5]<-"reflectance"
    
    ### The BLUP model is the same as it was for calculating BLUPs for grain yield.
    
    spectrum.blup.mod<-lmer(reflectance~(1|gid)+(1|trial)+(1|rep:trial)+(1|block:rep:trial), data=temp)
    
    ### The variance components can be extracted in the same way as with grain yield to 
    ### calculate broad-sense heritability.
    
    Vg<-data.frame(VarCorr(spectrum.blup.mod))$vcov[1]
    Ve<-data.frame(VarCorr(spectrum.blup.mod))$vcov[5]
    H2<-Vg/(Vg+(Ve/3))
    
    ### The broad-sense heritability can be put in the appropriate place within the 
    ### bands.H2 dataframe using which().
    
    bands.H2[which(bands.H2$band==bands[j] & bands.H2$date==names(spectra.list)[i]),"H2"]<-H2
    
    ### The BLUPs centered around the mean can also be calculated in the same way as with 
    ### grain yield and stored in a dataframe. 
    
    spectra.blups.temp<-ranef(spectrum.blup.mod)$gid+fixef(spectrum.blup.mod)
    spectra.blups.temp<-data.frame(rownames(spectra.blups.temp), spectra.blups.temp)
    colnames(spectra.blups.temp)<-c("gid", bands[j])
    
    ### The BLUPs can be merged based on the "gid" column with the existing dataframes 
    ### in the BLUPs list using the merge() function.

    spectra.blups.list[[i]]<-merge(spectra.blups.list[[i]], spectra.blups.temp, by="gid")

    ### Counters can be used to track progress
    
    print(i)
    print(j)
    
  }
  
}

View(spectra.blups.list)
head(bands.H2)

### To plot the broad-sense heritabilities for each band, the wavelengths need to be 
### extracted from the band labels as before. 

bands.H2$wavelength<-as.numeric(sapply(strsplit(as.character(bands.H2$band), "[.]"), "[", 2))
head(bands.H2)

### To plot different lines for each of the different dates, the color can be set to "date".

ggplot(bands.H2, aes(wavelength, H2, color=date)) + geom_line()

### To look at the spectral data for a single genotype, the which() function can be used
### within lapply() to pull out only observations for the gid 7305427. Using the "rbind"
### call within do.call() will stack these observations as rows in a data frame. The dates
### are listed as row names, so these can be added as a column using data.frame().

spectra.blups.sub<-do.call("rbind", lapply(spectra.blups.list, function(x) x[which(x$gid==7305427),]))
spectra.blups.sub<-data.frame(rownames(spectra.blups.sub), spectra.blups.sub)
colnames(spectra.blups.sub)[1]<-"date"

### The melt() function can be used to move the spectral bands from columns to rows.

spectra.blups.sub.melt<-melt(spectra.blups.sub)

### To plot the spectral data, the wavelengths need to be extracted from the band labels 
### as before. 

spectra.blups.sub.melt$wavelength<-as.numeric(sapply(strsplit(as.character(spectra.blups.sub.melt$variable), "[.]"), "[", 2))

### To plot different lines for each of the different dates, the color can be set to "date".

ggplot(spectra.blups.sub.melt, aes(wavelength, value, color=date)) + geom_line()

### To look at the hyperspectral data from multiple genotypes on a single data, the first 
### fifteen genotypes from the first date can be pulled out. As before, melt() moves the 
### spectral bands from columns to rows, and the wavelengths need to be extracted from the
### band labels. 

spectra.blups.sub.melt<-melt(spectra.blups.list[[1]][1:15,])
spectra.blups.sub.melt$wavelength<-as.numeric(sapply(strsplit(as.character(spectra.blups.sub.melt$variable), "[.]"), "[", 2))

### To plot different lines for each of the different genotypes, the color can be set to 
### "gid".

ggplot(spectra.blups.sub.melt, aes(wavelength, value, color=gid)) + geom_line()

### To prepare to calculate correlations between grain yield and each hyperspectral band
### on each date, a dataframe with a column for each band can be created. 

yield.corrs<-data.frame(bands)
colnames(yield.corrs)[1]<-"band"

for(i in 1:5){
  
  ### To calculate correlations, the first step is to merge the BLUPs for grain yield with
  ### the dataframe containing BLUPs for the hyperspectral bands from a single date. 
  
  spectra.yield.merged<-merge(yield.blups, spectra.blups.list[[i]], by="gid")
  
  ### The apply() function can be used to take the correlation between grain yield and each
  ### column of the dataframe, each containing a different band. 
  
  yield.corrs.temp<-apply(spectra.yield.merged[,-c(1:2)], 2, function(x) cor(spectra.yield.merged$grain.yield.blup, x))

  ### The correlations can be stored in a temporary dataframe with a column for the band 
  ### and another column for the correlation.
  
  yield.corrs.temp<-data.frame(names(yield.corrs.temp), yield.corrs.temp)
  
  ### The name of the column containing the correlations can be set to the date. 
  
  colnames(yield.corrs.temp)<-c("band", names(spectra.blups.list)[i])
  
  ### The correlations can be merged by band to the yield.corrs data frame previously 
  ### created. 
  
  yield.corrs<-merge(yield.corrs, yield.corrs.temp, by="band")
  
}

head(yield.corrs)

### To plot the correlations, melt() can be used to move them from columns to rows, and 
### the wavelengths can be extracted from the band labels.

yield.corrs.melt<-melt(yield.corrs)
colnames(yield.corrs.melt)
colnames(yield.corrs.melt)[c(2:3)]<-c("date", "yield.corr")

yield.corrs.melt$wavelength<-as.numeric(sapply(strsplit(as.character(yield.corrs.melt$band), "[.]"), "[", 2))
head(yield.corrs.melt)

### To plot different lines for each of the different dates, the color can be set to "date".

ggplot(yield.corrs.melt, aes(wavelength, yield.corr, color=date)) + geom_line()

### To prepare to create relationship matrices from the hyperspectral bands for each date, 
### an empty list with a length of five and the dates for names can be created as before.

spectra.rel.mat.list<-vector("list", 5)
names(spectra.rel.mat.list)<-names(spectra.blups.list)

### Relationship matrices can be created by taking the cross-product of the transposed 
### matrix of scaled BLUPs and dividing  by the number of spectral bands. The row and column
### names should be set to the gids. 

for(i in 1:length(spectra.rel.mat.list)){
  
  spectra.rel.mat.list[[i]]<-tcrossprod(scale(spectra.blups.list[[i]][,-1]))/ncol(spectra.blups.list[[i]][,-1])
  rownames(spectra.rel.mat.list[[i]])<-spectra.blups.list[[i]]$gid
  colnames(spectra.rel.mat.list[[i]])<-spectra.blups.list[[i]]$gid
  
}

spectra.rel.mat.list[[1]][1:10,1:10]

### To visualize part of one of the relationship matrices, the first 100 rows and columns 
### can be extracted. The melt() function will create a dataframe with two columns for 
### each pair of genotypes and a third column containing their covariance.

spectra.rel.mat.sub.melt<-melt(spectra.rel.mat.list[[1]][1:100,1:100])
spectra.rel.mat.sub.melt$Var1<-factor(spectra.rel.mat.sub.melt$Var1)
spectra.rel.mat.sub.melt$Var2<-factor(spectra.rel.mat.sub.melt$Var2)

### The relationship matrix can be visualized using geom_tile(). The "Var1" and "Var2" 
### columns created from the melt() function can be used as the x- and y-axes. The color
### of each tile can display the covariance by setting "fill" to the "value" column. To
### more easily distinguish differences, the "scale_fill_gradient2()" function can be used
### to specify which colors should be used for the low, mid, and high values. The last 
### line turns the gid labels along the x-axis on their side.

ggplot(spectra.rel.mat.sub.melt, aes(Var1, Var2, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2(low="red", mid="yellow", high="darkgreen") +
  theme(axis.text.x=element_text(angle=90))

### To prepare to test prediction models using each hyperspectral relationship matrix from
### each date, an empty list with a length of five and the dates for names can be created 
### as before.

pred.acc.list<-vector("list", 5)
names(pred.acc.list)<-names(spectra.list)

### The following nested loop will first go through the list of dataframes each containing
### hyperspectral relationship matrix from a single date.     

for(i in 1:length(spectra.rel.mat.list)){
  
  ### set.seed() is  called to make sure the different iterations use the same subsampling
  ### for each date.
  
  set.seed(90)
  
  ### Ten iterations of prediction are performed.
  
  for(j in 1:10){
    
    ### A random 20% of BLUEs for grain yield are set to NA within a temporary dataframe.
    
    temp<-yield.blues
    temp$grain.yield.blue[sample(1:nrow(temp), round(nrow(temp)/5))]<-NA
    
    ### kin.blup() from the "rrBLUP" package is then used to fit a prediction model using 
    ### the hyperspectral reflectance-derived relationship matrix.
    
    prediction.mod<-kin.blup(data=temp, geno="gid", pheno="grain.yield.blue", K=spectra.rel.mat.list[[i]])
    
    ### The predictions can be called using "$pred" with the model object. A dataframe
    ### with columns for "gid" and the predictions can be created. 
    
    preds<-prediction.mod$pred
    preds<-data.frame(names(preds), preds)
    colnames(preds)<-c("gid", "pred")
    
    ### This can be merged by gid to the "yield.blups" dataframe. To calculate the 
    ### prediction accuracy, only those gids for which the yield values were set to NA 
    ### should be considered. Those can be pulled out using which() to take out gids for 
    ### which the BLUE values are NA in the temporary dataframe. 
    
    preds.merged<-merge(yield.blups, preds, by="gid")
    preds.merged<-preds.merged[which(preds.merged$gid%in%temp[which(is.na(temp$grain.yield.blue)), "gid"]),]
    
    ### The prediction accuracy (Pearson's correlation between the prediction and the 
    ### observed BLUP for grain yield) can be appended to the appropriate list.
     
    pred.acc.list[[i]]<-append(pred.acc.list[[i]], cor(preds.merged$grain.yield.blup, preds.merged$pred))
    
    ### Counters can be printed to track progress. 
    
    print(i)
    print(j)
    
  }

}

pred.acc.list

### To prepare to plot the prediction accuracies, "rbind" can be called within do.call()
### to stack the prediction accuracies as rows in a single dataframe. melt() can then 
### be used such that each prediction accuracy appears in a single row. 

pred.acc<-melt(do.call("rbind", pred.acc.list))
head(pred.acc)

colnames(pred.acc)<-c("date", "iter", "pred.acc")

### geom_boxplot() can be used to visualize the prediction accuracies for the different
### dates. 

ggplot(pred.acc, aes(date, pred.acc)) + geom_boxplot()
