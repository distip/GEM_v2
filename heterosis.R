library(Matrix)
library(lme4)
library(reshape2)
library(rrBLUP)
library(ggplot2)
library(tidyverse)
library(zoom)
library(ggforce)
library(cowplot)

theme_set(theme_bw(10))

blups <- read.csv('spectra_blups.csv')
View(blups)

str(blups)

blups$Rep <- factor(blups$Rep)
blups$Block <- factor(blups$Block)
blups$year <- factor(blups$year)
blups$genotype <- factor(blups$genotype)
blups$note <- factor(blups$note)
blups$Trt <- factor(blups$Trt)
blups$ASD <- factor(blups$ASD)
blups$Group <- factor(blups$Group)
blups$rows <- factor(blups$rows)
blups$ranges <- factor(blups$ranges)
blups$PLOT.ID <- factor(blups$PLOT.ID)
blups$ASD  <- factor(blups$ASD)
blups$Calibration <- factor(blups$Calibration)

blups.melt <- melt(blups, id.vars = c('genotype', 'PLOT.ID', 'rows', 'ranges', 'Block', 'Rep', 'Trt', 'year', 'note', 'Calibration', 'ASD'))

View(blups.melt)
colnames(blups.melt)[12:13] <- c('band', 'reflectance')

blups.melt$band <- as.character(blups.melt$band)
blups.melt$band <- as.numeric(substr(blups.melt$band,2,5))

str(blups.melt)

blups.melt$note <-factor(blups.melt$note)

hybrid1<- blups.melt[which(blups.melt$genotype %in% c('BGEM-0109-N X B73', 'B73', 'BGEM-0109-N')),]
p1 <- ggplot(data = hybrid1, aes(x=band, y=reflectance, colour=genotype)) +
  geom_line() + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  facet_grid(cols  =vars(Trt))

hybrid2<- blups.melt[which(blups.melt$genotype %in% c('BGEM-0215-N X B73', 'B73', 'BGEM-0215-N')),]
p2 <- ggplot(data = hybrid2, aes(x=band, y=reflectance, colour=genotype)) +
  geom_line()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text.x = element_blank())+
  facet_grid(cols  =vars(Trt))


hybrid3<- blups.melt[which(blups.melt$genotype %in% c('BGEM-0223-N X B73', 'B73', 'BGEM-0223-N')),]
p3 <- ggplot(data = hybrid3, aes(x=band, y=reflectance, colour=genotype)) +
  geom_line()+
  theme(axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x=element_blank())+
  facet_grid(cols  =vars(Trt))


hybrid4<- blups.melt[which(blups.melt$genotype %in% c('BGEM-0134-S X Mo17', 'Mo17', 'BGEM-0134-S')),]
p4 <- ggplot(data = hybrid4, aes(x=band, y=reflectance, colour=genotype)) +
  geom_line()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank())+
  facet_grid(cols  =vars(Trt))

hybrid5<- blups.melt[which(blups.melt$genotype %in% c('BGEM-0264-S X Mo17', 'Mo17', 'BGEM-0264-S')),]
p5 <- ggplot(data = hybrid5, aes(x=band, y=reflectance, colour=genotype)) +
  geom_line()+
  theme(axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x=element_blank(),
        strip.text.x = element_blank())+
  facet_grid(cols  =vars(Trt))


hybrid6<- blups.melt[which(blups.melt$genotype %in% c('BGEM-0257-S X Mo17', 'Mo17', 'BGEM-0257-S')),]
p6 <- ggplot(data = hybrid6, aes(x=band, y=reflectance, colour=genotype)) +
  geom_line()+
  theme(strip.text.x = element_blank())+
  facet_grid(cols  =vars(Trt))

plot <- plot_grid(p3, p4, p5, p6, nrow = 4, ncol=1)
plot2<-grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob))

ggplot(data = hybrid[which(hybrid$Trt == 'LN'),], aes(x=band, y=reflectance, colour=genotype)) +
  geom_line()#+
  #facet_zoom(xlim = c(1380, 1410))

p2 <-ggplot(data = hybrid[hybrid$Trt == 'LN' & hybrid$band %in% c(1350:1450),], aes(x=band, y=reflectance, colour=genotype)) +
  geom_point(alpha=0.3)

p2


ggplot(data = hybrid[which(hybrid$Trt == 'LN'),], aes(x=band, y=reflectance, colour=genotype)) +
  stat_smooth(method='loess', span = 0.1, se=TRUE, aes(fill=genotype), alpha=0.5)




blups.mean <- blups.melt %>% group_by(note, band, Trt) %>% summarise(band=mean(band), reflectance=mean(reflectance))
View(blups.mean)


ggplot(data = blups.mean[which(blups.mean$Trt == 'HN'),]) +
  geom_line(aes(x = band, y = reflectance, color= note), alpha=0.4)

ggplot(data = blups.mean[which(blups.mean$Trt == 'LN'),],aes(x = band, y = reflectance, colour= note)) +
  stat_smooth(method='loess', span = 0.1, se=TRUE, aes(fill=note), alpha=0.3)



                    ########################   CALCULATION OF MID-PARENT HETEROSIS   #####################################
## heterosis calculation for low nitrogen condition
blups_LN <- blups[blups$Trt == 'LN' & blups$Rep == 1 ,]

het <- blups_LN %>% filter(grepl(' X ', genotype))

mid_LN <- expand.grid(bands= 350:2500, genotype = het$genotype)
mid_LN$heterosis <- NA

for(i in 1:length(het$genotype)){
  for(j in 13:length(colnames(het))){
    hybrid <- as.character(het$genotype[i])
    female <- strsplit(as.character(het$genotype[i]), ' X ')[[1]][1]
    male <- strsplit(as.character(het$genotype[i]), ' X ')[[1]][2]
    if(mean(blups_LN$genotype %in% c(female)) > 0 & mean(blups_LN$genotype %in% c(male) > 0 )) {
      female2 <- mean(blups_LN[which(blups_LN$genotype == female ), j])
      male2 <- mean(blups_LN[which(blups_LN$genotype == male ), j])
      hybrid2 <- mean(blups_LN[which(blups_LN$genotype == hybrid ), j])
      heterosis <- (hybrid2-mean(c(male2,female2)))/mean(c(male2,female2))*100
      mid_LN[mid_LN$bands == j+337 & mid_LN$genotype == hybrid , 'heterosis'] <- heterosis
      print(i)
      print(j)
    }
    else {
      print('pass')
    }
  }
 }

males <- c()
for(i in 1:length(mid_LN$genotype)){
  male <- strsplit(as.character(mid_LN$genotype[i]), ' X ')[[1]][2]
  males<- c(males,  male)
  print(i)
 }
mid_LN$male <- males
mid_LN <- na.omit(mid_LN)

p1 <-ggplot(data=mid_LN , aes(x=bands, y=heterosis, group=genotype, colour=male)) + 
  geom_line(size=0.4, alpha=0.6)+
  labs(title='Mid-Parent Heterosis in Low-Nitrogen', x='bands', y='heterosis %')
p1


p2 <- ggplot(data=mid_LN[mid_LN$male %in% c('B73', 'Mo17'),], aes(x=bands, y=heterosis, group=genotype, colour=male)) + 
  geom_line(size=0.4, alpha=0.6)+
  labs(title='Mid-Parent Heterosis in High-Nitrogen', x='bands', y='heterosis %')


data_mid_LN <- mid_LN[mid_LN$male %in% c('B73', 'Mo17'),] %>% group_by(male,bands) %>% 
  summarise(mean.heterosis = mean(heterosis, na.rm=TRUE), sd.heterosis = sd(heterosis, na.rm = TRUE), se.heterosis= sd(heterosis, na.rm=TRUE)/sqrt(length(heterosis)), 
            max = max(heterosis, na.rm = TRUE), min = min(heterosis, na.rm = TRUE))
View(data_mid_LN) 

p4 <- ggplot(data_mid_LN=data_mid_LN, aes(x=bands, group= male)) +
  geom_line(aes(y=mean.heterosis, color=male), size = 0.7)+
  geom_ribbon(aes(ymin=mean.heterosis-se.heterosis , ymax=mean.heterosis+se.heterosis , fill=male),alpha=0.4)+
  labs(title = 'Mid-parent Heterosis in LN', caption = '**Envelopes represent 1 se from the mean')
p4
######## heterosis calculation for high nitrogen condition ###########

## heterosis calculation for high nitrogen condition
blups_HN <- blups[blups$Trt == 'HN' & blups$Rep == 1 ,]

het <- blups_HN %>% filter(grepl(' X ', genotype))

mid.heterosis.val_HN <- expand.grid(bands= 350:2500, genotype = het$genotype)
mid.heterosis.val_HN$heterosis <- NA

for(i in 1:length(het$genotype)){
  for(j in 12:length(colnames(het))){
    hybrid <- as.character(het$genotype[i])
    female <- strsplit(as.character(het$genotype[i]), ' X ')[[1]][1]
    male <- strsplit(as.character(het$genotype[i]), ' X ')[[1]][2]
    if(mean(blups_HN$genotype %in% c(female)) > 0 & mean(blups_HN$genotype %in% c(male) > 0 )) {
      female2 <- mean(blups_HN[which(blups_HN$genotype == female ), j])
      male2 <- mean(blups_HN[which(blups_HN$genotype == male ), j])
      hybrid2 <- mean(blups_HN[which(blups_HN$genotype == hybrid ), j])
      heterosis <- (hybrid2-mean(c(male2,female2)))/mean(c(male2,female2))*100
      mid.heterosis.val_HN[mid.heterosis.val_HN$bands == j+338 & mid.heterosis.val_HN$genotype == hybrid , 'heterosis'] <- heterosis
      print(i)
      print(j)
    }
    else
      print('pass')
  }
}

males <- c()
for(i in 1:length(mid.heterosis.val_HN$genotype)){
  male <- strsplit(as.character(mid.heterosis.val_HN$genotype[i]), ' X ')[[1]][2]
  males<- c(males,  male)
  print(i)
}
mid.heterosis.val_HN$male <- males
mid.heterosis.val_HN <- na.omit(mid.heterosis.val_HN)
mid.heterosis.val_HN$male <- factor(mid.heterosis.val_HN$male)
copy <- mid.heterosis.val_HN

p2 <- ggplot(data=mid.heterosis.val_HN[mid.heterosis.val_HN$male %in% c('B73', 'Mo17'),], aes(x=bands, y=heterosis, group=genotype, colour=male)) + 
  geom_line(size=0.4, alpha=0.6)+
  labs(title='Mid-Parent Heterosis in High-Nitrogen', x='bands', y='heterosis %')

p3 <- ggplot(data=mid.heterosis.val_HN[mid.heterosis.val_HN$male %in% c('B73', 'Mo17'),], aes(x=bands, y=heterosis, colour=male)) +
  stat_smooth(method='loess', span = 0.1, se=FALSE, aes(fill=male), alpha=0.5)

dat2 <- mid.heterosis.val_HN[mid.heterosis.val_HN$male %in% c('B73', 'Mo17'),] %>% group_by(male,bands) %>% 
  summarise(mean.heterosis = mean(heterosis, na.rm=TRUE), sd.heterosis = sd(heterosis, na.rm = TRUE), se.heterosis= sd(heterosis, na.rm=TRUE)/sqrt(length(heterosis)), 
            max = max(heterosis, na.rm = TRUE), min = min(heterosis, na.rm = TRUE))


p4 <- ggplot(data=data, aes(x=bands, group= male)) +
  geom_line(aes(y=mean.heterosis, color=male), size = 0.7)+
  geom_ribbon(aes(ymin=mean.heterosis-sd.heterosis , ymax=mean.heterosis+sd.heterosis , fill=male),alpha=0.4)+
  labs(title = 'Mid-parent Heterosis in HN', caption = '**Envelopes represent 1 sd from the mean')
p4


### merging 2 datafrmae to make a common plot

data$Trt <- c(rep('LN', length(data$bands)))
data_lower_HN$Trt <- c(rep('HN', length(data_lower_HN$bands)))

plot_data <- data %>% full_join(data2)

p6 <- ggplot(data=plot_data, aes(x=bands, group= male)) +
  geom_line(data= plot_data[plot_data$Trt == 'HN',], aes( x= bands,y=mean.heterosis, linetype=male, color= Trt),size = 0.7)+
  geom_ribbon(data= plot_data[plot_data$Trt == 'HN',], aes(ymin=mean.heterosis-se.heterosis , ymax=mean.heterosis+se.heterosis),alpha=0.2)+
  geom_line(data= plot_data[plot_data$Trt == 'LN',], aes( x= bands,y=mean.heterosis,linetype=male, color = Trt), size = 0.7)+
  geom_ribbon(data= plot_data[plot_data$Trt == 'LN',], aes(ymin=mean.heterosis-se.heterosis , ymax=mean.heterosis+se.heterosis),alpha=0.2)+
  labs(title = 'Mid-parent Heterosis under HN and LN conditions', caption = '**Envelopes represent 1 se from the mean')

p6







                ########################   CALCULATION OF lower-PARENT HETEROSIS   #####################################
## heterosis calculation for low nitrogen condition
blups_LN <- blups[blups$Trt == 'LN' & blups$Rep == 1 ,]

het <- blups_LN %>% filter(grepl(' X ', genotype))

heterosis.val_LN <- expand.grid(bands= 350:2500, genotype = het$genotype)
heterosis.val_LN$heterosis <- NA

for(i in 1:length(het$genotype)){
  for(j in 13:length(colnames(het))){
    hybrid <- as.character(het$genotype[i])
    female <- strsplit(as.character(het$genotype[i]), ' X ')[[1]][1]
    male <- strsplit(as.character(het$genotype[i]), ' X ')[[1]][2]
    if(mean(blups_LN$genotype %in% c(female)) > 0 & mean(blups_LN$genotype %in% c(male) > 0 )) {
      female2 <- mean(blups_LN[which(blups_LN$genotype == female ), j])
      male2 <- mean(blups_LN[which(blups_LN$genotype == male ), j])
      lower <- min(c(female2, male2))
      hybrid2 <- mean(blups_LN[which(blups_LN$genotype == hybrid ), j])
      heterosis <- (lower-hybrid2)/lower*100
      heterosis.val_LN[heterosis.val_LN$bands == j+337 & heterosis.val_LN$genotype == hybrid , 'heterosis'] <- heterosis
      print(i)
      print(j)
    }
    else
      print('pass')
  }
}

males <- c()
for(i in 1:length(heterosis.val_LN$genotype)){
  male <- strsplit(as.character(heterosis.val_LN$genotype[i]), ' X ')[[1]][2]
  males<- c(males,  male)
  print(i)
}
heterosis.val_LN$male <- males
heterosis.val_LN <- na.omit(heterosis.val_LN)

p1 <-ggplot(data=heterosis.val_LN , aes(x=bands, y=heterosis, group=genotype, colour=male)) + 
  geom_line(size=0.4, alpha=0.6)+
  labs(title='Lower-Parent Heterosis in Low-Nitrogen', x='bands', y='heterosis %')
p1


p2 <- ggplot(data=heterosis.val_LN[heterosis.val_LN$male %in% c('B73', 'Mo17'),], aes(x=bands, y=heterosis, group=genotype, colour=male)) + 
  geom_line(size=0.4, alpha=0.6)+
  labs(title='Lower-Parent Heterosis in Low-Nitrogen', x='bands', y='heterosis %')


data <- heterosis.val_LN[heterosis.val_LN$male %in% c('B73', 'Mo17'),] %>% group_by(male,bands) %>% 
  summarise(mean.heterosis = mean(heterosis, na.rm=TRUE), sd.heterosis = sd(heterosis, na.rm = TRUE), se.heterosis= sd(heterosis, na.rm=TRUE)/sqrt(length(heterosis)), 
            max = max(heterosis, na.rm = TRUE), min = min(heterosis, na.rm = TRUE))
View(data) 

p4 <- ggplot(data=data, aes(x=bands, group= male)) +
  geom_line(aes(y=mean.heterosis, color=male), size = 0.7)+
  geom_ribbon(aes(ymin=mean.heterosis-se.heterosis , ymax=mean.heterosis+se.heterosis , fill=male),alpha=0.4)+
  labs(title = 'Lower-parent Heterosis in LN', caption = '**Envelopes represent 1 sd from the mean')
p4
######## lower parent heterosis calculation for high nitrogen condition ###########

## heterosis calculation for high nitrogen condition
blups_HN <- blups[blups$Trt == 'HN' & blups$Rep == 1 ,]

het <- blups_HN %>% filter(grepl(' X ', genotype))

lower_HN <- expand.grid(bands= 350:2500, genotype = het$genotype)
lower_HN$heterosis <- NA

for(i in 1:length(het$genotype)){
  for(j in 13:length(colnames(het))){
    hybrid <- as.character(het$genotype[i])
    female <- strsplit(as.character(het$genotype[i]), ' X ')[[1]][1]
    male <- strsplit(as.character(het$genotype[i]), ' X ')[[1]][2]
    if(mean(blups_HN$genotype %in% c(female)) > 0 & mean(blups_HN$genotype %in% c(male) > 0 )) {
      female2 <- mean(blups_HN[which(blups_HN$genotype == female ), j])
      male2 <- mean(blups_HN[which(blups_HN$genotype == male ), j])
      lower <- min(c(female2, male2))
      hybrid2 <- mean(blups_HN[which(blups_HN$genotype == hybrid ), j])
      heterosis <- (lower-hybrid2)/lower*100
      lower_HN[lower_HN$bands == j+337 & lower_HN$genotype == hybrid , 'heterosis'] <- heterosis
      print(i)
      print(j)
    }
    else
      print('pass')
  }
}

males <- c()
for(i in 1:length(lower_HN$genotype)){
  male <- strsplit(as.character(lower_HN$genotype[i]), ' X ')[[1]][2]
  males<- c(males,  male)
  print(i)
}
lower_HN$male <- males
lower_HN <- na.omit(lower_HN)
lower_HN$male <- factor(lower_HN$male)
copy <- lower_HN


p2 <- ggplot(data=lower_HN[lower_HN$male %in% c('B73', 'Mo17'),], aes(x=bands, y=heterosis, group=genotype, colour=male)) + 
  geom_line(size=0.4, alpha=0.6)+
  labs(title='Lower-Parent Heterosis in High-Nitrogen', x='bands', y='heterosis %')
p2

p3 <- ggplot(data=lower_HN[lower_HN$male %in% c('B73', 'Mo17'),], aes(x=bands, y=heterosis, colour=male)) +
  stat_smooth(method='loess', span = 0.1, se=FALSE, aes(fill=male), alpha=0.5)

data_lower_HN <- lower_HN[lower_HN$male %in% c('B73', 'Mo17'),] %>% group_by(male,bands) %>% 
  summarise(mean.heterosis = mean(heterosis, na.rm=TRUE), sd.heterosis = sd(heterosis, na.rm = TRUE), se.heterosis= sd(heterosis, na.rm=TRUE)/sqrt(length(heterosis)), 
            max = max(heterosis, na.rm = TRUE), min = min(heterosis, na.rm = TRUE))


p5 <- ggplot(data=data_lower_HN, aes(x=bands, group= male)) +
  geom_line(aes(y=mean.heterosis, color=male), size = 0.7)+
  geom_ribbon(aes(ymin=mean.heterosis-se.heterosis , ymax=mean.heterosis+se.heterosis , fill=male),alpha=0.4)+
  labs(title = 'Lower-parent Heterosis in HN', caption = '**Envelopes represent 1 se from the mean')
p5



### merging 2 datafrmae to make a common plot

data$Trt <- c(rep('LN', length(data$bands)))
data2$Trt <- c(rep('HN', length(data$bands)))

plot_data <- data %>% full_join(data2)

p6 <- ggplot(data=plot_data, aes(x=bands, group= male)) +
  geom_line(data= plot_data[plot_data$Trt == 'HN',], aes( x= bands,y=mean.heterosis, linetype=male, color= Trt),size = 0.7)+
  geom_ribbon(data= plot_data[plot_data$Trt == 'HN',], aes(ymin=mean.heterosis-se.heterosis , ymax=mean.heterosis+se.heterosis),alpha=0.2)+
  geom_line(data= plot_data[plot_data$Trt == 'LN',], aes( x= bands,y=mean.heterosis,linetype=male, color = Trt), size = 0.7)+
  geom_ribbon(data= plot_data[plot_data$Trt == 'LN',], aes(ymin=mean.heterosis-se.heterosis , ymax=mean.heterosis+se.heterosis),alpha=0.2)+
  labs(title = 'Lower-parent Heterosis under HN and LN conditions', caption = '**Envelopes represent 1 se from the mean')

p6



########################   CALCULATION OF Better-PARENT HETEROSIS   #####################################
## heterosis calculation for low nitrogen condition
blups_LN <- blups[blups$Trt == 'LN' & blups$Rep == 1 ,]

het <- blups_LN %>% filter(grepl(' X ', genotype))

heterosis.val_LN <- expand.grid(bands= 350:2500, genotype = het$genotype)
heterosis.val_LN$heterosis <- NA

for(i in 1:length(het$genotype)){
  for(j in 13:length(colnames(het))){
    hybrid <- as.character(het$genotype[i])
    female <- strsplit(as.character(het$genotype[i]), ' X ')[[1]][1]
    male <- strsplit(as.character(het$genotype[i]), ' X ')[[1]][2]
    if(mean(blups_LN$genotype %in% c(female)) > 0 & mean(blups_LN$genotype %in% c(male) > 0 )) {
      female2 <- mean(blups_LN[which(blups_LN$genotype == female ), j])
      male2 <- mean(blups_LN[which(blups_LN$genotype == male ), j])
      better <- max(c(female2, male2))
      hybrid2 <- mean(blups_LN[which(blups_LN$genotype == hybrid ), j])
      heterosis <- (hybrid2-better)/better*100
      heterosis.val_LN[heterosis.val_LN$bands == j+337 & heterosis.val_LN$genotype == hybrid , 'heterosis'] <- heterosis
      print(i)
      print(j)
    }
    else
      print('pass')
  }
}

males <- c()
for(i in 1:length(heterosis.val_LN$genotype)){
  male <- strsplit(as.character(heterosis.val_LN$genotype[i]), ' X ')[[1]][2]
  males<- c(males,  male)
  print(i)
}
heterosis.val_LN$male <- males
heterosis.val_LN <- na.omit(heterosis.val_LN)

p1 <-ggplot(data=heterosis.val_LN , aes(x=bands, y=heterosis, group=genotype, colour=male)) + 
  geom_line(size=0.4, alpha=0.6)+
  labs(title='Better-Parent Heterosis in Low-Nitrogen', x='bands', y='heterosis %')
p1


p2 <- ggplot(data=heterosis.val_LN[heterosis.val_LN$male %in% c('B73', 'Mo17'),], aes(x=bands, y=heterosis, group=genotype, colour=male)) + 
  geom_line(size=0.4, alpha=0.6)+
  labs(title='Lower-Parent Heterosis in Low-Nitrogen', x='bands', y='heterosis %')


data <- heterosis.val_LN[heterosis.val_LN$male %in% c('B73', 'Mo17'),] %>% group_by(male,bands) %>% 
  summarise(mean.heterosis = mean(heterosis, na.rm=TRUE), sd.heterosis = sd(heterosis, na.rm = TRUE), se.heterosis= sd(heterosis, na.rm=TRUE)/sqrt(length(heterosis)), 
            max = max(heterosis, na.rm = TRUE), min = min(heterosis, na.rm = TRUE))
View(data) 

p4 <- ggplot(data=data, aes(x=bands, group= male)) +
  geom_line(aes(y=mean.heterosis, color=male), size = 0.7)+
  geom_ribbon(aes(ymin=mean.heterosis-se.heterosis , ymax=mean.heterosis+se.heterosis , fill=male),alpha=0.4)+
  labs(title = 'Better-parent Heterosis in LN', caption = '**Envelopes represent 1 se from the mean')
p4
######## heterosis calculation for high nitrogen condition ###########

## heterosis calculation for high nitrogen condition
blups_HN <- blups[blups$Trt == 'HN' & blups$Rep == 1 ,]

het <- blups_HN %>% filter(grepl(' X ', genotype))

mid.heterosis.val_HN <- expand.grid(bands= 350:2500, genotype = het$genotype)
mid.heterosis.val_HN$heterosis <- NA

for(i in 1:length(het$genotype)){
  for(j in 13:length(colnames(het))){
    hybrid <- as.character(het$genotype[i])
    female <- strsplit(as.character(het$genotype[i]), ' X ')[[1]][1]
    male <- strsplit(as.character(het$genotype[i]), ' X ')[[1]][2]
    if(mean(blups_HN$genotype %in% c(female)) > 0 & mean(blups_HN$genotype %in% c(male) > 0 )) {
      female2 <- mean(blups_HN[which(blups_HN$genotype == female ), j])
      male2 <- mean(blups_HN[which(blups_HN$genotype == male ), j])
      better <- max(c(female2, male2))
      hybrid2 <- mean(blups_HN[which(blups_HN$genotype == hybrid ), j])
      heterosis <- (hybrid2-better)/better*100
      mid.heterosis.val_HN[mid.heterosis.val_HN$bands == j+337 & mid.heterosis.val_HN$genotype == hybrid , 'heterosis'] <- heterosis
      print(i)
      print(j)
    }
    else
      print('pass')
  }
}

males <- c()
for(i in 1:length(mid.heterosis.val_HN$genotype)){
  male <- strsplit(as.character(mid.heterosis.val_HN$genotype[i]), ' X ')[[1]][2]
  males<- c(males,  male)
  print(i)
}
mid.heterosis.val_HN$male <- males
mid.heterosis.val_HN <- na.omit(mid.heterosis.val_HN)
mid.heterosis.val_HN$male <- factor(mid.heterosis.val_HN$male)
copy <- mid.heterosis.val_HN


p2 <- ggplot(data=mid.heterosis.val_HN[mid.heterosis.val_HN$male %in% c('B73', 'Mo17'),], aes(x=bands, y=heterosis, group=genotype, colour=male)) + 
  geom_line(size=0.4, alpha=0.6)+
  labs(title='Lower-Parent Heterosis in High-Nitrogen', x='bands', y='heterosis %')
p2

p3 <- ggplot(data=mid.heterosis.val_HN[mid.heterosis.val_HN$male %in% c('B73', 'Mo17'),], aes(x=bands, y=heterosis, colour=male)) +
  stat_smooth(method='loess', span = 0.1, se=FALSE, aes(fill=male), alpha=0.5)

data2 <- mid.heterosis.val_HN[mid.heterosis.val_HN$male %in% c('B73', 'Mo17'),] %>% group_by(male,bands) %>% 
  summarise(mean.heterosis = mean(heterosis, na.rm=TRUE), sd.heterosis = sd(heterosis, na.rm = TRUE), se.heterosis= sd(heterosis, na.rm=TRUE)/sqrt(length(heterosis)), 
            max = max(heterosis, na.rm = TRUE), min = min(heterosis, na.rm = TRUE))


p5 <- ggplot(data=data2, aes(x=bands, group= male)) +
  geom_line(aes(y=mean.heterosis, color=male), size = 0.7)+
  geom_ribbon(aes(ymin=mean.heterosis-se.heterosis , ymax=mean.heterosis+se.heterosis , fill=male),alpha=0.4)+
  labs(title = 'Better-parent Heterosis in HN', caption = '**Envelopes represent 1 se from the mean')
p5



### merging 2 datafrmae to make a common plot

data$Trt <- c(rep('LN', length(data$bands)))
data2$Trt <- c(rep('HN', length(data$bands)))

plot_data <- data %>% full_join(data2)

p6 <- ggplot(data=plot_data, aes(x=bands, group= male)) +
  geom_line(data= plot_data[plot_data$Trt == 'HN',], aes( x= bands,y=mean.heterosis, linetype=male, color= Trt),size = 0.7)+
  geom_ribbon(data= plot_data[plot_data$Trt == 'HN',], aes(ymin=mean.heterosis-se.heterosis , ymax=mean.heterosis+se.heterosis),alpha=0.2)+
  geom_line(data= plot_data[plot_data$Trt == 'LN',], aes( x= bands,y=mean.heterosis,linetype=male, color = Trt), size = 0.7)+
  geom_ribbon(data= plot_data[plot_data$Trt == 'LN',], aes(ymin=mean.heterosis-se.heterosis , ymax=mean.heterosis+se.heterosis),alpha=0.2)+
  labs(title = 'Better-parent Heterosis under HN and LN conditions', caption = '**Envelopes represent 1 se from the mean')

p6
