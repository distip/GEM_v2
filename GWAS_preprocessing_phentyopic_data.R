library(tidyverse)

GWAS_data <- read.csv('blups_mergedv2_with_gwas_ids.csv')

View(GWAS_data)

GWAS_data_HN <- GWAS_data[GWAS_data$Trt == 'HN', ]

GWAS <- subset(GWAS_data_HN, select = -c(X, genotype, new_GID, PLOT.ID, rows, ranges,  Block, Rep, Trt, year, note, Group, Calibration,  ASD))

write.csv(GWAS, './pheno_HN.csv', row.names = FALSE)

View(GWAS)



                                      ######################### Data preperation for CHL  #####################

library(tidyverse)

GWAS_data <- read.csv('blues_CHL_with_gwas_ids_only_inbreds.csv')

View(GWAS_data)

GWAS_data_HN <- GWAS_data

GWAS <- subset(GWAS_data_HN, select = c(UID ,CHL))

write.csv(GWAS, './pheno_CHL_blues_only_inbreds.csv', row.names = FALSE)

View(GWAS)
