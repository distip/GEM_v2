
data <- read.csv('Raw_spectrum_merged_predicted_CHL.csv')

CHL <- data %>% select('PLOT.ID', 'CHL')

View(CHL)

data2 <- read.csv('Raw_spectrum_merged.csv')

first_cols <- data2 %>% select('PLOT.ID' ,'genotype','Trt', 'Rep' , 'Group')


merged <- merge(first_cols, CHL , all.x = TRUE)

View(merged)

write.csv(merged, './Raw_CHL_predicted', row.names = FALSE)
