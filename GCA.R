library(tidyverse)
library(gpbStat)
library(agricolae)

blups <- read.csv('spectra_blups.csv')
hybrids <- blups[blups$note == 'Hybrid',]

str(hybrids)


hybrids <- hybrids %>% filter(grepl(' X ', genotype))
View(hybrids)

female <- c()
male <- c()
for (i in 1:length(hybrids$genotype)){
  h <- as.character(hybrids$genotype[i])
  f <- strsplit(as.character(hybrids$genotype[i]), ' X ')[[1]][1]
  m <- strsplit(as.character(hybrids$genotype[i]), ' X ')[[1]][2]
  female <- c(female, f)
  male <- c(male, m)
  
}

hybrids <- hybrids %>% mutate(female = female , male = male, .before = note)
grouped <- hybrids %>% group_by(genotype, Rep) %>% summarise(trait = mean(X350), male= male, female =female)
grouped <- grouped %>% group_by(genotype, )

output_LN <- with(grouped, lineXtester(Rep, female, male, trait))

