library(Matrix)
library(lme4)
library(reshape2)
library(rrBLUP)
library(ggplot2)
library(tidyverse)
library(viridisLite)
library(viridis)
library(grid)

data <- read_csv('blups_merged_v2')
data <- data.frame(data)

data$Rep <- factor(data$Rep)
data$Block <- factor(data$Block)
data$year <- factor(data$year)
data$genotype <- factor(data$genotype)
data$note <- factor(data$note)
data$Trt <- factor(data$Trt)
data$ASD <- factor(data$ASD)
data$Group <- factor(data$Group)
data$rows <- factor(data$rows)
data$ranges <- factor(data$ranges)
data$PLOT.ID <- factor(data$PLOT.ID)
data$ASD  <- factor(data$ASD)
data$Calibration <- factor(data$Calibration)

View(data)


ggplot(data= data %>% select(genotype,X550,Trt))+
  geom_point(aes(x=Trt, y=X550))+
  geom_line(aes(x=Trt, y=X550, group=genotype, color='red'))

                                ####### plasticity calculation ################

bands <- (colnames(spectra)[-c(1:12)])
bands_df <-as.data.frame(colnames(spectra)[-c(1:12)])
bands

spectral.plasticity.list <- data.frame(levels(data$genotype))
colnames(spectral.plasticity.list) <- 'genotype'

for(j in 1:length(bands)){
  for(i in 1:length(spectral.plasticity.list$genotype)){
    ### The hyperspectral data from just one N can be stored in a temporary variable.
    temp <- data
    geno <- spectral.plasticity.list$genotype[i]
    geno_hn <- unique(data[data$genotype == geno & data$Trt == 'HN', bands[j]])
    geno_ln <- unique(data[data$genotype == geno & data$Trt == 'LN', bands[j]])
    if(length(geno_hn) != 0  & length(geno_ln) != 0 ){
      plasticity <- (geno_ln - geno_hn) / geno_hn
      spectral.plasticity.list[spectral.plasticity.list$genotype == geno, bands[j]] <- plasticity
    }
    else{
      print('cokerler')
    }
    print(j)
  }
}

bview(spectral.plasticity.list)  


data_group <- data %>% select(genotype, Group)
data_group <- data_group[!duplicated(data_group$genotype), ]

merged <- merge(data_group, spectral.plasticity.list)

View(merged)

merged_melt <- melt(data= merged, id.vars = c('genotype', 'Group'))             

colnames(merged_melt)[c(3,4)] <- c('band', 'plasticity')

merged_melt$band <- as.character(merged_melt$band)
merged_melt$band<- as.numeric(substr(merged_melt$band,2,5))

p <- ggplot(data=merged_melt, aes(x=band, y=plasticity, color=Group))+ #[merged_melt$Group == 'Inbred', ]
  geom_point(size=0.5, alpha=0.04)+
  guides(colour=guide_legend(override.aes=list(alpha= 1)))
p + theme_bw()

View(spectral.plasticity.list)

spectra_plasticity <- spectral.plasticity.list %>% mutate_if(is.numeric, round, digits=3)
View(spectra_plasticity)
merged <- merge(data_group, spectra_plasticity)

write.csv(merged, 'spectra_plasticity.csv', row.names = FALSE)