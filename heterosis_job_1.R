library(lme4)
library(reshape2)
library(rrBLUP)
library(ggplot2)
library(tidyverse)
library(zoom)
library(ggforce)
library(cowplot)


blups <- read.csv("spectra_comb_blups_v2.csv")

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




########################   CALCULATION OF MID-PARENT HETEROSIS   #####################################
## heterosis calculation for low nitrogen condition
blups_LN <- blups[blups$Trt == "LN" & blups$Rep == 2 ,]

het <- blups_LN %>% filter(grepl(" X ", genotype))

mid_LN <- expand.grid(bands= 350:2500, genotype = het$genotype)
mid_LN$heterosis <- NA

for(i in 1:length(het$genotype)){
  for(j in 14:length(colnames(het))){
    hybrid <- as.character(het$genotype[i])
    female <- strsplit(as.character(het$genotype[i]), " X ")[[1]][1]
    male <- strsplit(as.character(het$genotype[i]), " X ")[[1]][2]
    if(mean(blups_LN$genotype %in% c(female)) > 0 & mean(blups_LN$genotype %in% c(male) > 0 )) {
      female2 <- mean(blups_LN[which(blups_LN$genotype == female ), j])
      male2 <- mean(blups_LN[which(blups_LN$genotype == male ), j])
      hybrid2 <- mean(blups_LN[which(blups_LN$genotype == hybrid ), j])
      heterosis <- (hybrid2-mean(c(male2,female2)))/mean(c(male2,female2))*100
      mid_LN[mid_LN$bands == j+336 & mid_LN$genotype == hybrid , "heterosis"] <- heterosis
      print(i)
      print(j)
    }
    else {
      print("pass")
    }
  }
}

males <- c()
for(i in 1:length(mid_LN$genotype)){
  male <- strsplit(as.character(mid_LN$genotype[i]), " X ")[[1]][2]
  males<- c(males,  male)
  print(i)
}
mid_LN$male <- males
mid_LN <- na.omit(mid_LN)

write_csv(mid_LN, 'mid_LN_v2.csv')

data_mid_LN <- mid_LN[mid_LN$male %in% c("B73", "Mo17"),] %>% group_by(male,bands) %>% 
  summarise(mean.heterosis = mean(heterosis, na.rm=TRUE), sd.heterosis = sd(heterosis, na.rm = TRUE), se.heterosis= sd(heterosis, na.rm=TRUE)/sqrt(length(heterosis)), 
            max = max(heterosis, na.rm = TRUE), min = min(heterosis, na.rm = TRUE))


######## mid parent heterosis calculation for high nitrogen condition ###########

## heterosis calculation for high nitrogen condition
blups_HN <- blups[blups$Trt == "HN" & blups$Rep == 2 ,]

het <- blups_HN %>% filter(grepl(" X ", genotype))

mid_HN <- expand.grid(bands= 350:2500, genotype = het$genotype)
mid_HN$heterosis <- NA

for(i in 1:length(het$genotype)){
  for(j in 13:length(colnames(het))){
    hybrid <- as.character(het$genotype[i])
    female <- strsplit(as.character(het$genotype[i]), " X ")[[1]][1]
    male <- strsplit(as.character(het$genotype[i]), " X ")[[1]][2]
    if(mean(blups_HN$genotype %in% c(female)) > 0 & mean(blups_HN$genotype %in% c(male) > 0 )) {
      female2 <- mean(blups_HN[which(blups_HN$genotype == female ), j])
      male2 <- mean(blups_HN[which(blups_HN$genotype == male ), j])
      hybrid2 <- mean(blups_HN[which(blups_HN$genotype == hybrid ), j])
      heterosis <- (hybrid2-mean(c(male2,female2)))/mean(c(male2,female2))*100
      mid_HN[mid_HN$bands == j+337 & mid_HN$genotype == hybrid , "heterosis"] <- heterosis
      print(i)
      print(j)
    }
    else
      print("pass")
  }
}

males <- c()
for(i in 1:length(mid_HN$genotype)){
  male <- strsplit(as.character(mid_HN$genotype[i]), " X ")[[1]][2]
  males<- c(males,  male)
  print(i)
}
mid_HN$male <- males
mid_HN <- na.omit(mid_HN)

write_csv(mid_HN, 'mid_HN.csv')
data_mid_HN <- mid_HN[mid_HN$male %in% c("B73", "Mo17"),] %>% group_by(male,bands) %>% 
  summarise(mean.heterosis = mean(heterosis, na.rm=TRUE), sd.heterosis = sd(heterosis, na.rm = TRUE), se.heterosis= sd(heterosis, na.rm=TRUE)/sqrt(length(heterosis)), 
            max = max(heterosis, na.rm = TRUE), min = min(heterosis, na.rm = TRUE))

########################   CALCULATION OF lower-PARENT HETEROSIS   #####################################
## heterosis calculation for low nitrogen condition
blups_LN <- blups[blups$Trt == "LN" & blups$Rep == 1 ,]

het <- blups_LN %>% filter(grepl(" X ", genotype))

low_LN <- expand.grid(bands= 350:2500, genotype = het$genotype)
low_LN$heterosis <- NA

for(i in 1:length(het$genotype)){
  for(j in 13:length(colnames(het))){
    hybrid <- as.character(het$genotype[i])
    female <- strsplit(as.character(het$genotype[i]), " X ")[[1]][1]
    male <- strsplit(as.character(het$genotype[i]), " X ")[[1]][2]
    if(mean(blups_LN$genotype %in% c(female)) > 0 & mean(blups_LN$genotype %in% c(male) > 0 )) {
      female2 <- mean(blups_LN[which(blups_LN$genotype == female ), j])
      male2 <- mean(blups_LN[which(blups_LN$genotype == male ), j])
      lower <- min(c(female2, male2))
      hybrid2 <- mean(blups_LN[which(blups_LN$genotype == hybrid ), j])
      heterosis <- (lower-hybrid2)/lower*100
      low_LN[low_LN$bands == j+337 & low_LN$genotype == hybrid , "heterosis"] <- heterosis
      print(i)
      print(j)
    }
    else
      print("pass")
  }
}

males <- c()
for(i in 1:length(low_LN$genotype)){
  male <- strsplit(as.character(low_LN$genotype[i]), " X ")[[1]][2]
  males<- c(males,  male)
  print(i)
}
low_LN$male <- males
low_LN <- na.omit(low_LN)

write_csv(low_LN, 'low_LN.csv')

data_low_LN <- low_LN[low_LN$male %in% c("B73", "Mo17"),] %>% group_by(male,bands) %>% 
  summarise(mean.heterosis = mean(heterosis, na.rm=TRUE), sd.heterosis = sd(heterosis, na.rm = TRUE), se.heterosis= sd(heterosis, na.rm=TRUE)/sqrt(length(heterosis)), 
            max = max(heterosis, na.rm = TRUE), min = min(heterosis, na.rm = TRUE))

######## lower parent heterosis calculation for high nitrogen condition ###########

## heterosis calculation for high nitrogen condition
blups_HN <- blups[blups$Trt == "HN" & blups$Rep == 1 ,]

het <- blups_HN %>% filter(grepl(" X ", genotype))

low_HN <- expand.grid(bands= 350:2500, genotype = het$genotype)
low_HN$heterosis <- NA

for(i in 1:length(het$genotype)){
  for(j in 13:length(colnames(het))){
    hybrid <- as.character(het$genotype[i])
    female <- strsplit(as.character(het$genotype[i]), " X ")[[1]][1]
    male <- strsplit(as.character(het$genotype[i]), " X ")[[1]][2]
    if(mean(blups_HN$genotype %in% c(female)) > 0 & mean(blups_HN$genotype %in% c(male) > 0 )) {
      female2 <- mean(blups_HN[which(blups_HN$genotype == female ), j])
      male2 <- mean(blups_HN[which(blups_HN$genotype == male ), j])
      lower <- min(c(female2, male2))
      hybrid2 <- mean(blups_HN[which(blups_HN$genotype == hybrid ), j])
      heterosis <- (lower-hybrid2)/lower*100
      low_HN[low_HN$bands == j+337 & low_HN$genotype == hybrid , "heterosis"] <- heterosis
      print(i)
      print(j)
    }
    else
      print("pass")
  }
}

males <- c()
for(i in 1:length(low_HN$genotype)){
  male <- strsplit(as.character(low_HN$genotype[i]), " X ")[[1]][2]
  males<- c(males,  male)
  print(i)
}
low_HN$male <- males
low_HN <- na.omit(low_HN)

write_csv(low_HN, 'low_HN.csv')

data_low_HN <- low_HN[low_HN$male %in% c("B73", "Mo17"),] %>% group_by(male,bands) %>% 
  summarise(mean.heterosis = mean(heterosis, na.rm=TRUE), sd.heterosis = sd(heterosis, na.rm = TRUE), se.heterosis= sd(heterosis, na.rm=TRUE)/sqrt(length(heterosis)), 
            max = max(heterosis, na.rm = TRUE), min = min(heterosis, na.rm = TRUE))



########################   CALCULATION OF Better-PARENT HETEROSIS   #####################################
## heterosis calculation for low nitrogen condition
blups_LN <- blups[blups$Trt == "LN" & blups$Rep == 1 ,]

het <- blups_LN %>% filter(grepl(" X ", genotype))

better_LN <- expand.grid(bands= 350:2500, genotype = het$genotype)
better_LN$heterosis <- NA

for(i in 1:length(het$genotype)){
  for(j in 13:length(colnames(het))){
    hybrid <- as.character(het$genotype[i])
    female <- strsplit(as.character(het$genotype[i]), " X ")[[1]][1]
    male <- strsplit(as.character(het$genotype[i]), " X ")[[1]][2]
    if(mean(blups_LN$genotype %in% c(female)) > 0 & mean(blups_LN$genotype %in% c(male) > 0 )) {
      female2 <- mean(blups_LN[which(blups_LN$genotype == female ), j])
      male2 <- mean(blups_LN[which(blups_LN$genotype == male ), j])
      better <- max(c(female2, male2))
      hybrid2 <- mean(blups_LN[which(blups_LN$genotype == hybrid ), j])
      heterosis <- (hybrid2-better)/better*100
      better_LN[better_LN$bands == j+337 & better_LN$genotype == hybrid , "heterosis"] <- heterosis
      print(i)
      print(j)
    }
    else
      print("pass")
  }
}

males <- c()
for(i in 1:length(better_LN$genotype)){
  male <- strsplit(as.character(better_LN$genotype[i]), " X ")[[1]][2]
  males<- c(males,  male)
  print(i)
}
better_LN$male <- males
better_LN <- na.omit(better_LN)

write_csv(better_LN, 'better_LN.csv')

data_better_LN <- better_LN[better_LN$male %in% c("B73", "Mo17"),] %>% group_by(male,bands) %>% 
  summarise(mean.heterosis = mean(heterosis, na.rm=TRUE), sd.heterosis = sd(heterosis, na.rm = TRUE), se.heterosis= sd(heterosis, na.rm=TRUE)/sqrt(length(heterosis)), 
            max = max(heterosis, na.rm = TRUE), min = min(heterosis, na.rm = TRUE))

######## heterosis calculation for high nitrogen condition ###########

## heterosis calculation for high nitrogen condition
blups_HN <- blups[blups$Trt == "HN" & blups$Rep == 1 ,]

het <- blups_HN %>% filter(grepl(" X ", genotype))

better_HN <- expand.grid(bands= 350:2500, genotype = het$genotype)
better_HN$heterosis <- NA

for(i in 1:length(het$genotype)){
  for(j in 13:length(colnames(het))){
    hybrid <- as.character(het$genotype[i])
    female <- strsplit(as.character(het$genotype[i]), " X ")[[1]][1]
    male <- strsplit(as.character(het$genotype[i]), " X ")[[1]][2]
    if(mean(blups_HN$genotype %in% c(female)) > 0 & mean(blups_HN$genotype %in% c(male) > 0 )) {
      female2 <- mean(blups_HN[which(blups_HN$genotype == female ), j])
      male2 <- mean(blups_HN[which(blups_HN$genotype == male ), j])
      better <- max(c(female2, male2))
      hybrid2 <- mean(blups_HN[which(blups_HN$genotype == hybrid ), j])
      heterosis <- (hybrid2-better)/better*100
      better_HN[better_HN$bands == j+337 & better_HN$genotype == hybrid , "heterosis"] <- heterosis
      print(i)
      print(j)
    }
    else
      print("pass")
  }
}

males <- c()
for(i in 1:length(better_HN$genotype)){
  male <- strsplit(as.character(better_HN$genotype[i]), " X ")[[1]][2]
  males<- c(males,  male)
  print(i)
}
better_HN$male <- males
better_HN <- na.omit(better_HN)

write_csv(better_HN, 'better_HN.csv')

data_better_HN <- better_HN[better_HN$male %in% c("B73", "Mo17"),] %>% group_by(male,bands) %>% 
  summarise(mean.heterosis = mean(heterosis, na.rm=TRUE), sd.heterosis = sd(heterosis, na.rm = TRUE), se.heterosis= sd(heterosis, na.rm=TRUE)/sqrt(length(heterosis)), 
            max = max(heterosis, na.rm = TRUE), min = min(heterosis, na.rm = TRUE))


                                    ############################### MERGING ALL DATAFRAMES ###########################


data_mid_HN$Trt <- c(rep('HN', length(data_mid_HN$bands)))
data_mid_HN$Htype <- c(rep('mid', length(data_mid_HN$bands)))

data_mid_LN$Trt <- c(rep('LN', length(data_mid_LN$bands)))
data_mid_LN$Htype <- c(rep('mid', length(data_mid_LN$bands)))

data_better_HN$Trt <- c(rep('HN', length(data_better_HN$bands)))
data_better_HN$Htype <- c(rep('better', length(data_better_HN$bands)))

data_better_LN$Trt <- c(rep('LN', length(data_better_LN$bands)))
data_better_LN$Htype <- c(rep('better', length(data_better_LN$bands)))

data_low_HN$Trt <- c(rep('HN', length(data_low_HN$bands)))
data_low_HN$Htype <- c(rep('low', length(data_low_HN$bands)))

data_low_LN$Trt <- c(rep('LN', length(data_low_LN$bands)))
data_low_LN$Htype <- c(rep('low', length(data_low_LN$bands)))



plot_data <- data_mid_HN %>% full_join(data_mid_LN) %>% full_join(data_better_HN) %>% full_join(data_better_LN) %>% full_join(data_low_HN) %>%
  full_join(data_low_LN)

write.csv(plot_data, 'heterosis_data.csv')

sub_plot_data = plot_data[plot_data$Htype == 'mid',]

theme_set(theme_bw(10))
p6 <- ggplot(data=sub_plot_data, aes(x=bands, group= male)) +
  geom_line(data= sub_plot_data[sub_plot_data$Trt == 'HN',], aes( x= bands,y=mean.heterosis, linetype=male, color= Trt),size = 0.7)+
  geom_ribbon(data= sub_plot_data[sub_plot_data$Trt == 'HN',], aes(ymin=mean.heterosis-se.heterosis , ymax=mean.heterosis+se.heterosis),alpha=0.2)+
  geom_line(data= sub_plot_data[sub_plot_data$Trt == 'LN',], aes( x= bands,y=mean.heterosis,linetype=male, color = Trt), size = 0.7)+
  geom_ribbon(data= sub_plot_data[sub_plot_data$Trt == 'LN',], aes(ymin=mean.heterosis-se.heterosis , ymax=mean.heterosis+se.heterosis),alpha=0.2)+
  scale_x_continuous(breaks =c(350, 550, 710, 1445,1890, 2020, 2410, 2500))+
  theme(axis.text.x = element_text(angle=45))+
  labs(title = 'Mid-parent Heterosis under HN and LN conditions', caption = '**Envelopes represent 1 se from the mean')

p6 +
  annotate('rect', xmin=350, xmax=750, ymin=-10, ymax=10, alpha= 0.1, fill='blue')+
  annotate('rect', xmin=750, xmax=2500, ymin=-10, ymax=10, alpha= 0.1, fill='yellow')+
  annotate('text', x=550, y=9, label= 'Environmentally \n Driven')+
  annotate('text', x=1500, y=9, label = 'Genetically Driven')+
  geom_vline(xintercept = 550, size=0.2)+
  geom_vline(xintercept = 710, size=0.2)+
  geom_vline(xintercept = 1445, size=0.2)+
  geom_vline(xintercept = 1890, size=0.2)+
  geom_vline(xintercept = 2020, size=0.2)+
  geom_vline(xintercept = 2410, size=0.2)


analysis_mid_LN <- mid_LN[mid_LN$bands %in%  c('550', '710', '1445', '1890', '2020', '2410'),] # %>% distinct(genotype, .keep_all = TRUE)
analysis_mid_HN <- mid_HN[mid_HN$bands %in% c('550', '710', '1445', '1890', '2020', '2410') ,] # %>% distinct(genotype, .keep_all = TRUE)
analysis_mid_HN$Trt <- c(rep('HN', length(analysis_mid_HN$bands)))
analysis_mid_LN$Trt <- c(rep('LN', length(analysis_mid_LN$bands)))

analysis_mid_merge <- analysis_mid_HN %>% full_join(analysis_mid_LN)
analysis_mid_merge$bands <- factor(analysis_mid_merge$bands)

box_mid <- ggplot(data = analysis_mid_merge, aes(x = bands, y = heterosis, fill= male)) + 
  geom_boxplot(width=0.1, position = position_dodgenudge(width=1))+
  geom_boxplot(width=0.1, position= position_nudge(x = -0.2),  aes(fill= Trt))+
  #annotate('text' , x=0.8 , y= 22, label='p-value < 2.2e-16')+
  labs(title = 'MPH values of Spectra in 6 hotspots')
  

box_mid

t.test(analysis_mid_HN[analysis_mid_HN$bands == '1890', 'heterosis'], analysis_mid_LN[analysis_mid_LN$bands == '1890', 'heterosis'], var.equal = TRUE)





