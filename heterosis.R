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

heterosis.val_LN <- expand.grid(bands= 350:2500, genotype = het$genotype)
heterosis.val_LN$heterosis <- NA

for(i in 1:length(het$genotype)){
  for(j in 12:length(colnames(het))){
    hybrid <- as.character(het$genotype[i])
    female <- strsplit(as.character(het$genotype[i]), ' X ')[[1]][1]
    male <- strsplit(as.character(het$genotype[i]), ' X ')[[1]][2]
    if(mean(blups_LN$genotype %in% c(female)) > 0 & mean(blups_LN$genotype %in% c(male) > 0 )) {
      female2 <- mean(blups_LN[which(blups_LN$genotype == female ), j])
      male2 <- mean(blups_LN[which(blups_LN$genotype == male ), j])
      hybrid2 <- mean(blups_LN[which(blups_LN$genotype == hybrid ), j])
      heterosis <- (hybrid2-mean(c(male2,female2)))/mean(c(male2,female2))*100
      heterosis.val_LN[heterosis.val_LN$bands == j+338 & heterosis.val_LN$genotype == hybrid , 'heterosis'] <- heterosis
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
  labs(title='Mid-Parent Heterosis in Low-Nitrogen', x='bands', y='heterosis %')
p1


p2 <- ggplot(data=heterosis.val_LN[heterosis.val_LN$male %in% c('B73', 'Mo17'),], aes(x=bands, y=heterosis, group=genotype, colour=male)) + 
  geom_line(size=0.4, alpha=0.6)+
  labs(title='Mid-Parent Heterosis in High-Nitrogen', x='bands', y='heterosis %')


data <- heterosis.val_LN[heterosis.val_LN$male %in% c('B73', 'Mo17'),] %>% group_by(male,bands) %>% 
  summarise(mean.heterosis = mean(heterosis, na.rm=TRUE), sd.heterosis = sd(heterosis, na.rm = TRUE), se.heterosis= sd(heterosis, na.rm=TRUE)/sqrt(length(heterosis)), 
            max = max(heterosis, na.rm = TRUE), min = min(heterosis, na.rm = TRUE))
View(data) 

p4 <- ggplot(data=data, aes(x=bands, group= male)) +
  geom_line(aes(y=mean.heterosis, color=male), size = 0.7)+
  geom_ribbon(aes(ymin=mean.heterosis-sd.heterosis , ymax=mean.heterosis+sd.heterosis , fill=male),alpha=0.4)+
  labs(title = 'Mid-parent Heterosis in LN', caption = '**Envelopes represent 1 sd from the mean')
p4
######## heterosis calculation for high nitrogen condition ###########

## heterosis calculation for high nitrogen condition
blups_HN <- blups[blups$Trt == 'HN' & blups$Rep == 1 ,]

het <- blups_HN %>% filter(grepl(' X ', genotype))

heterosis.val_HN <- expand.grid(bands= 350:2500, genotype = het$genotype)
heterosis.val_HN$heterosis <- NA

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
      heterosis.val_HN[heterosis.val_HN$bands == j+338 & heterosis.val_HN$genotype == hybrid , 'heterosis'] <- heterosis
      print(i)
      print(j)
    }
    else
      print('pass')
  }
}

males <- c()
for(i in 1:length(heterosis.val_HN$genotype)){
  male <- strsplit(as.character(heterosis.val_HN$genotype[i]), ' X ')[[1]][2]
  males<- c(males,  male)
  print(i)
}
heterosis.val_HN$male <- males
heterosis.val_HN <- na.omit(heterosis.val_HN)
heterosis.val_HN$male <- factor(heterosis.val_HN$male)
copy <- heterosis.val_HN

p2 <- ggplot(data=heterosis.val_HN[heterosis.val_HN$male %in% c('B73', 'Mo17'),], aes(x=bands, y=heterosis, group=genotype, colour=male)) + 
  geom_line(size=0.4, alpha=0.6)+
  labs(title='Mid-Parent Heterosis in High-Nitrogen', x='bands', y='heterosis %')

p3 <- ggplot(data=heterosis.val_HN[heterosis.val_HN$male %in% c('B73', 'Mo17'),], aes(x=bands, y=heterosis, colour=male)) +
  stat_smooth(method='loess', span = 0.1, se=FALSE, aes(fill=male), alpha=0.5)

data <- heterosis.val_HN[heterosis.val_HN$male %in% c('B73', 'Mo17'),] %>% group_by(male,bands) %>% 
  summarise(mean.heterosis = mean(heterosis, na.rm=TRUE), sd.heterosis = sd(heterosis, na.rm = TRUE), se.heterosis= sd(heterosis, na.rm=TRUE)/sqrt(length(heterosis)), 
            max = max(heterosis, na.rm = TRUE), min = min(heterosis, na.rm = TRUE))


p4 <- ggplot(data=data, aes(x=bands, group= male)) +
  geom_line(aes(y=mean.heterosis, color=male), size = 0.7)+
  geom_ribbon(aes(ymin=mean.heterosis-sd.heterosis , ymax=mean.heterosis+sd.heterosis , fill=male),alpha=0.4)+
  labs(title = 'Mid-parent Heterosis in HN', caption = '**Envelopes represent 1 sd from the mean')
p4

model <- lm(heterosis ~ 1 , heterosis.val_HN)
confint(model,level = 0.95)







                ########################   CALCULATION OF Better-PARENT HETEROSIS   #####################################
## heterosis calculation for low nitrogen condition
blups_LN <- blups[blups$Trt == 'LN' & blups$Rep == 1 ,]

het <- blups_LN %>% filter(grepl(' X ', genotype))

heterosis.val_LN <- expand.grid(bands= 350:2500, genotype = het$genotype)
heterosis.val_LN$heterosis <- NA

for(i in 1:length(het$genotype)){
  for(j in 12:length(colnames(het))){
    hybrid <- as.character(het$genotype[i])
    female <- strsplit(as.character(het$genotype[i]), ' X ')[[1]][1]
    male <- strsplit(as.character(het$genotype[i]), ' X ')[[1]][2]
    if(mean(blups_LN$genotype %in% c(female)) > 0 & mean(blups_LN$genotype %in% c(male) > 0 )) {
      female2 <- mean(blups_LN[which(blups_LN$genotype == female ), j])
      male2 <- mean(blups_LN[which(blups_LN$genotype == male ), j])
      better <- max(c(female2, male2))
      hybrid2 <- mean(blups_LN[which(blups_LN$genotype == hybrid ), j])
      heterosis <- (hybrid2-better)/better*100
      heterosis.val_LN[heterosis.val_LN$bands == j+338 & heterosis.val_LN$genotype == hybrid , 'heterosis'] <- heterosis
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
  labs(title='Better-Parent Heterosis in Low-Nitrogen', x='bands', y='heterosis %')


data <- heterosis.val_LN[heterosis.val_LN$male %in% c('B73', 'Mo17'),] %>% group_by(male,bands) %>% 
  summarise(mean.heterosis = mean(heterosis, na.rm=TRUE), sd.heterosis = sd(heterosis, na.rm = TRUE), se.heterosis= sd(heterosis, na.rm=TRUE)/sqrt(length(heterosis)), 
            max = max(heterosis, na.rm = TRUE), min = min(heterosis, na.rm = TRUE))
View(data) 

p4 <- ggplot(data=data, aes(x=bands, group= male)) +
  geom_line(aes(y=mean.heterosis, color=male), size = 0.7)+
  geom_ribbon(aes(ymin=mean.heterosis-sd.heterosis , ymax=mean.heterosis+sd.heterosis , fill=male),alpha=0.4)+
  labs(title = 'Better-parent Heterosis in LN', caption = '**Envelopes represent 1 sd from the mean')
p4
######## heterosis calculation for high nitrogen condition ###########

## heterosis calculation for high nitrogen condition
blups_HN <- blups[blups$Trt == 'HN' & blups$Rep == 1 ,]

het <- blups_HN %>% filter(grepl(' X ', genotype))

heterosis.val_HN <- expand.grid(bands= 350:2500, genotype = het$genotype)
heterosis.val_HN$heterosis <- NA

for(i in 1:length(het$genotype)){
  for(j in 12:length(colnames(het))){
    hybrid <- as.character(het$genotype[i])
    female <- strsplit(as.character(het$genotype[i]), ' X ')[[1]][1]
    male <- strsplit(as.character(het$genotype[i]), ' X ')[[1]][2]
    if(mean(blups_HN$genotype %in% c(female)) > 0 & mean(blups_HN$genotype %in% c(male) > 0 )) {
      female2 <- mean(blups_HN[which(blups_HN$genotype == female ), j])
      male2 <- mean(blups_HN[which(blups_HN$genotype == male ), j])
      better <- max(c(female2, male2))
      hybrid2 <- mean(blups_HN[which(blups_HN$genotype == hybrid ), j])
      heterosis <- (hybrid2-better)/better*100
      heterosis.val_HN[heterosis.val_HN$bands == j+338 & heterosis.val_HN$genotype == hybrid , 'heterosis'] <- heterosis
      print(i)
      print(j)
    }
    else
      print('pass')
  }
}

males <- c()
for(i in 1:length(heterosis.val_HN$genotype)){
  male <- strsplit(as.character(heterosis.val_HN$genotype[i]), ' X ')[[1]][2]
  males<- c(males,  male)
  print(i)
}
heterosis.val_HN$male <- males
heterosis.val_HN <- na.omit(heterosis.val_HN)
heterosis.val_HN$male <- factor(heterosis.val_HN$male)
copy <- heterosis.val_HN

p2 <- ggplot(data=heterosis.val_HN[heterosis.val_HN$male %in% c('B73', 'Mo17'),], aes(x=bands, y=heterosis, group=genotype, colour=male)) + 
  geom_line(size=0.4, alpha=0.6)+
  labs(title='Better-Parent Heterosis in High-Nitrogen', x='bands', y='heterosis %')
p2

p3 <- ggplot(data=heterosis.val_HN[heterosis.val_HN$male %in% c('B73', 'Mo17'),], aes(x=bands, y=heterosis, colour=male)) +
  stat_smooth(method='loess', span = 0.1, se=FALSE, aes(fill=male), alpha=0.5)

data <- heterosis.val_HN[heterosis.val_HN$male %in% c('B73', 'Mo17'),] %>% group_by(male,bands) %>% 
  summarise(mean.heterosis = mean(heterosis, na.rm=TRUE), sd.heterosis = sd(heterosis, na.rm = TRUE), se.heterosis= sd(heterosis, na.rm=TRUE)/sqrt(length(heterosis)), 
            max = max(heterosis, na.rm = TRUE), min = min(heterosis, na.rm = TRUE))


p4 <- ggplot(data=data, aes(x=bands, group= male)) +
  geom_line(aes(y=mean.heterosis, color=male), size = 0.7)+
  geom_ribbon(aes(ymin=mean.heterosis-sd.heterosis , ymax=mean.heterosis+sd.heterosis , fill=male),alpha=0.4)+
  labs(title = 'Better-parent Heterosis in HN', caption = '**Envelopes represent 1 sd from the mean')
p4

model <- lm(heterosis ~ 1 , heterosis.val_HN)
confint(model,level = 0.95)