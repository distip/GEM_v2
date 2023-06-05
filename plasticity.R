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

View(spectral.plasticity.list)  


data_group <- data %>% select(genotype, Group)
data_group <- data_group[!duplicated(data_group$genotype), ]

merged <- merge(data_group, spectral.plasticity.list)

View(merged)

merged_melt <- melt(data= merged, id.vars = c('genotype', 'Group'))             

colnames(merged_melt)[c(3,4)] <- c('band', 'plasticity')

merged_melt$band <- as.character(merged_melt$band)
merged_melt$band<- as.numeric(substr(merged_melt$band,2,5))

plasticity_plot <- merged_melt %>% group_by(Group, band) %>% summarise(mean.pl = mean(plasticity, na.rm=TRUE), sd.pl = sd(plasticity, na.rm=TRUE))


p2 <- ggplot(data=plasticity_plot, aes(x=band)) +
  geom_line(aes(x=band, y=mean.pl))+
  geom_ribbon(aes(ymin=mean.pl-sd.pl , ymax= mean.pl + sd.pl),alpha=0.2)+
  facet_wrap(~ Group)+
  labs(y='Plasticity', x='', caption = '**Envelopes represent 1 sd from the mean')+
  theme_bw(12)
p2

p <- ggplot(data=merged_melt, 
            aes(x=band, y=plasticity, color=Group))+ #[merged_melt$Group == 'Inbred', ]
  geom_point(size=0.5, alpha=0.1)+
  guides(colour=guide_legend(override.aes=list(alpha= 1)))+
  facet_wrap(~ Group)
p + theme_bw(16)+
  geom_vline(xintercept = 560, size=0.2, linetype= 'dashed')+
  geom_vline(xintercept = 710, size=0.2, linetype= 'dashed')+
  geom_vline(xintercept = 1445, size=0.2, linetype= 'dashed')+
  geom_vline(xintercept = 1890, size=0.2, linetype= 'dashed')+
  geom_vline(xintercept = 2010, size=0.2, linetype= 'dashed')+
  geom_vline(xintercept = 2420, size=0.2, linetype= 'dashed')+
  scale_x_continuous(breaks =c(350, 560, 710, 1445,1890, 2010, 2420, 2500))+
  theme(axis.text.x = element_text(angle=45))


write.csv(merged, 'spectra_plasticity_from_blups_v1.csv', row.names = FALSE)



                #### Plotting Plasticity for biotraits ####
data <- read.csv('blups_all_predicted.csv')

data_HN <- data %>% select(genotype, SLA_HN, LWC_HN, X..N_HN, X..P_HN, X..K_HN, X..S_HN,
                           X..Ca_HN, X..Mg_HN, ppm.Zn_HN, ppm.Fe_HN, ppm.Mn_HN,
                           ppm.Cu_HN, ppm.B_HN, ppm.Mo_HN, SLA_HN, LWC_HN)  
#data_HN <- mutate(data_HN , Trt=rep('HN', times= length(data_HN$genotype)))
colnames(data_HN) <- gsub('_HN', '', colnames(data_HN))
data_HN <- data_HN %>% melt(id.vars= c('genotype'))
colnames(data_HN)[2:3] <- c('traits', 'HN' )

data_LN <- data %>% select(genotype, SLA_LN, LWC_LN, X..N_LN, X..P_LN, X..K_LN, X..S_LN,
                           X..Ca_LN, X..Mg_LN, ppm.Zn_LN, ppm.Fe_LN, ppm.Mn_LN,
                           ppm.Cu_LN, ppm.B_LN, ppm.Mo_LN, SLA_LN, LWC_LN)  
#data_LN <- mutate(data_LN , Trt=rep('LN', times= length(data_LN$genotype)))
colnames(data_LN) <- gsub('_LN', '', colnames(data_LN))
data_LN <- data_LN %>% melt(id.vars= c('genotype'))
colnames(data_LN)[2:3] <- c('traits', 'LN') 

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
  facet_wrap('traits', scales = 'free')+
  labs(title = 'Nitrogen Responses of Biochemical Traits')+
  ylab('Traits')+
  xlab('')


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

