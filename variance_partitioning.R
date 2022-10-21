library(Matrix)
library(lme4)
library(reshape2)
library(rrBLUP)
library(ggplot2)
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
spectra$Group <- factor(spectra$Group)
spectra$rows <- factor(spectra$rows)
spectra$ranges <- factor(spectra$ranges)
spectra$PLOT.ID <- factor(spectra$PLOT.ID)
spectra$ASD  <- factor(spectra$ASD)
spectra$Calibration <- factor(spectra$Calibration)

spectra <- subset(spectra, select = -c(X, Unnamed..0))

spectra.hybrids <- spectra[spectra$Group == 'Inbred',]


bands <- (colnames(spectra.hybrids)[-c(1:13)])
bands_df <-as.data.frame(colnames(spectra.hybrids)[-c(1:13)])
bands 
bands_df$H2 <- NA
colnames(bands_df) <- c('band', 'H2')
bands.H2 <- bands_df

spectra.hybridsl.blups.list <- data.frame(levels(spectra.hybrids$genotype))

colnames(spectra.hybridsl.blups.list) <- 'genotype'



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
  temp <- spectra.hybrids
  temp<-temp[, which(colnames(spectra.hybrids) %in% c('PLOT.ID','genotype', 'Block', 'Trt', 'Rep', 'ASD', bands[i]))]
  colnames(temp)[7] <- 'reflectance'
  m <- lmer(reflectance ~Trt + (Trt|genotype), data=temp)
  extractVarsLmm(m)
  LN <- extractVarsLmm(m)[[1]]/ sum(extractVarsLmm(m))
  HN <- extractVarsLmm(m)[[2]]/ sum(extractVarsLmm(m))
  Nitrate <- extractVarsLmm(m)[[3]]/ sum(extractVarsLmm(m))
  Plasticity <- extractVarsLmm(m)[[4]]/ sum(extractVarsLmm(m))
  Res <- extractVarsLmm(m)[[5]]/ sum(extractVarsLmm(m))
  values <-  c(LN, HN, Nitrate, Plasticity, Res)
  var.part.list[[i]] <- values
  print(i)
}

var.part.list.melt <- melt(var.part.list)
var.part.list.melt$source <- rep(c('LN', 'HN', 'Nitrate', 'Plasticity', 'Residual'), 10750/5)
var.part.list.melt$source <- factor(var.part.list.melt$source)
colnames(var.part.list.melt) <- c('values', 'band', 'source')
var.part.list.melt$band <- as.numeric(var.part.list.melt$band)




ggplot(var.part.list.melt, aes(fill= source, y=values, x=band)) +
  geom_bar(position='fill', stat='identity', width = 1, alpha=0.8)+
  theme_classic()+
  theme(legend.background = element_rect(fill = '#FFCC66',color='grey50',  size=1))+
  scale_fill_brewer(palette= 'Set3')+
  labs(title = 'Variance Partitioning of Leaf Spectrum' , subtitle = 'lmer(reflectance ~ Trt + genotype + (Trt|genotype), data=temp)' , y ='%', x='wavelengths')
  








