library(tidyverse)

GWAS_data <- read.csv('BGEM_extra_pheno_for_GWAS_blues.csv')

View(GWAS_data)

GWAS_data <- GWAS_data[GWAS_data$Trt == 'LN', ]

GWAS <- subset(GWAS_data, select = c(UID, leaf_length, leaf_width, ear_height, flag_leaf, plant_height))

View(GWAS)

write.csv(GWAS, './pheno_LN_blues.csv', row.names = FALSE)

View(GWAS)



                                      ######################### Data preperation for CHL  #####################

library(tidyverse)

GWAS_data <- read.csv('blues_CHL_with_gwas_ids_only_inbreds.csv')

View(GWAS_data)

GWAS_data_HN <- GWAS_data

GWAS <- subset(GWAS_data_HN, select = c(UID ,CHL))

write.csv(GWAS, './pheno_CHL_blues_only_inbreds.csv', row.names = FALSE)

View(GWAS)
