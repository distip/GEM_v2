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
                #### Plotting Plasticity for ind wawelengths ####
data_HN <- data[data$Trt == 'HN',] %>% select(genotype, Group, X560, X710, X1445, X1890, X2010, X2420) %>%
  melt(id.vars= c('genotype', 'Group'))  
colnames(data_HN)[3:4] <- c('band', 'HN')
data_LN <- data[data$Trt == 'LN',] %>% select(genotype, Group, X560, X710, X1445, X1890, X2010, X2420) %>%
  melt(id.vars= c('genotype', 'Group'))
colnames(data_LN)[3:4] <- c('band', 'LN')
plot_data <- merge(data_HN, data_LN)
plot_data <- mutate(plot_data, plasticity = (LN-HN)/HN,
                    plasticity = case_when(plasticity < -0.5 ~ -0.5,
                                           plasticity > 0.5 ~ 0.5,
                                           TRUE ~ plasticity)) %>% 
  pivot_longer(c(HN, LN), names_to = 'Trt', values_to='reflectance') 

plot_data %>% 
  ggplot(aes(x=Trt, y=reflectance))+
  geom_line(aes(group=genotype, color=plasticity), show.legend = TRUE, size=0.6)+
  scale_color_gradient2(name='Plasticity',
                        high='#FF0000',
                        mid ='#C0C0C0',
                        low = '#0000FF',
                        limits= c(-0.5, 0.5),
                        breaks = c(0.5, 0.25, 0, -0.25, -0.5),
                        labels = c('>0.50', '0.25', '0' , '-0.25', '<-0.50'))+
   geom_boxplot(aes(x=Trt, y=reflectance), width=0.2, outlier.shape = NA, show.legend = FALSE)+
   geom_point(aes(x=Trt, y=reflectance), size=0.3)+
   theme_bw(16)+
   facet_wrap('band', scales = 'free')+
   labs(title = 'Nitrogen Responses of Selected Wavelengths')+
   ylab('Reflectance')+
   xlab('')



                                ####### plasticity calculation ################
spectra <- data
bands <- (colnames(spectra)[-c(1:13)])
bands_df <-as.data.frame(colnames(spectra)[-c(1:13)])
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

view(spectral.plasticity.list)  


data_group <- data %>% select(genotype, Group)
data_group <- data_group[!duplicated(data_group$genotype), ]

merged <- merge(data_group, spectral.plasticity.list)

View(merged)

merged_melt <- melt(data= merged, id.vars = c('genotype', 'Group'))             

colnames(merged_melt)[c(3,4)] <- c('band', 'plasticity')

merged_melt$band <- as.character(merged_melt$band)
merged_melt$band<- as.numeric(substr(merged_melt$band,2,5))

p <- ggplot(data=merged_melt, 
            aes(x=band, y=plasticity, color=Group))+ #[merged_melt$Group == 'Inbred', ]
  geom_point(size=0.5, alpha=0.1)+
  guides(colour=guide_legend(override.aes=list(alpha= 1)))
p + theme_bw(16)+
  geom_vline(xintercept = 560, size=0.2)+
  geom_vline(xintercept = 710, size=0.2)+
  geom_vline(xintercept = 1445, size=0.2)+
  geom_vline(xintercept = 1890, size=0.2)+
  geom_vline(xintercept = 2010, size=0.2)+
  geom_vline(xintercept = 2420, size=0.2)+
  scale_x_continuous(breaks =c(350, 560, 710, 1445,1890, 2010, 2420, 2500))+
  theme(axis.text.x = element_text(angle=45))


write.csv(merged, 'spectra_plasticity_from_blups_v1.csv', row.names = FALSE)
