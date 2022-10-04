library(tidyverse)
library(gpbStat)

hybrids <- blups[blups$note == 'Hybrid',]

str(hybrids)


hybrids <- hybrids %>% filter(grepl(' X ', genotype))
View(hybrids)

female <- c()
male <- c()
for (i in 1:length(hybrids$note)){
  h <- as.character(hybrids$genotype[i])
  f <- strsplit(as.character(het$genotype[i]), ' X ')[[1]][1]
  m <- strsplit(as.character(het$genotype[i]), ' X ')[[1]][2]
  female <- c(female, f)
  male <- c(male, m)
  
}
male

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