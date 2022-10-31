library(lme4)
library(reshape2)
library(rrBLUP)
library(ggplot2)
library(tidyverse)
library(zoom)
library(ggforce)
library(cowplot)


blups <- read.csv("spectra_blups.csv")

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
blups_LN <- blups[blups$Trt == "LN" & blups$Rep == 1 ,]

het <- blups_LN %>% filter(grepl(" X ", genotype))

mid_LN <- expand.grid(bands= 350:2500, genotype = het$genotype)
mid_LN$heterosis <- NA

for(i in 1:length(het$genotype)){
  for(j in 13:length(colnames(het))){
    hybrid <- as.character(het$genotype[i])
    female <- strsplit(as.character(het$genotype[i]), " X ")[[1]][1]
    male <- strsplit(as.character(het$genotype[i]), " X ")[[1]][2]
    if(mean(blups_LN$genotype %in% c(female)) > 0 & mean(blups_LN$genotype %in% c(male) > 0 )) {
      female2 <- mean(blups_LN[which(blups_LN$genotype == female ), j])
      male2 <- mean(blups_LN[which(blups_LN$genotype == male ), j])
      hybrid2 <- mean(blups_LN[which(blups_LN$genotype == hybrid ), j])
      heterosis <- (hybrid2-mean(c(male2,female2)))/mean(c(male2,female2))*100
      mid_LN[mid_LN$bands == j+337 & mid_LN$genotype == hybrid , "heterosis"] <- heterosis
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


data_mid_LN <- mid_LN[mid_LN$male %in% c("B73", "Mo17"),] %>% group_by(male,bands) %>% 
  summarise(mean.heterosis = mean(heterosis, na.rm=TRUE), sd.heterosis = sd(heterosis, na.rm = TRUE), se.heterosis= sd(heterosis, na.rm=TRUE)/sqrt(length(heterosis)), 
            max = max(heterosis, na.rm = TRUE), min = min(heterosis, na.rm = TRUE))

######## mid parent heterosis calculation for high nitrogen condition ###########

## heterosis calculation for high nitrogen condition
blups_HN <- blups[blups$Trt == "HN" & blups$Rep == 1 ,]

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


data_better_HN <- better_HN[better_HN$male %in% c("B73", "Mo17"),] %>% group_by(male,bands) %>% 
  summarise(mean.heterosis = mean(heterosis, na.rm=TRUE), sd.heterosis = sd(heterosis, na.rm = TRUE), se.heterosis= sd(heterosis, na.rm=TRUE)/sqrt(length(heterosis)), 
            max = max(heterosis, na.rm = TRUE), min = min(heterosis, na.rm = TRUE))
