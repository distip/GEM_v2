library(FactoMineR)
library(ggplot2)
library(factoextra)
library(tidyverse)


spectra <- read.csv("Raw_spectrum_merged")
View(spectra)

str(spectra)

spectra$Rep <- factor(spectra$Rep)
spectra$Block <- factor(spectra$Block)
spectra$year <- factor(spectra$year)
spectra$genotype <- factor(spectra$genotype)
spectra$note <- factor(spectra$note)
spectra$Trt <- factor(spectra$Trt)
spectra$ASD <- factor(spectra$ASD)


blups <- read.csv('spectra_blups.csv')

str(blups)

blups$genotype <- as.factor(blups$genotype)
blups$PLOT.ID <- as.factor(blups$PLOT.ID)
blups$rows <- as.factor(blups$rows)
blups$ranges <- as.factor(blups$ranges)
blups$Block <- as.factor(blups$Block)
blups$Rep <- as.factor(blups$Rep)
blups$Trt <- as.factor(blups$Trt)
blups$year <- as.factor(blups$year)
blups$note <-as.factor(blups$note)
blups$Calibration <- as.factor(blups$Calibration)
blups$ASD <- as.factor(blups$ASD)

## adding 'hybrid or inbred' column
group <- c()
for(i in 1:length(blups$note)){
  note <- blups$note[i]
  if(note == 'Hybrid'){
    group <- c(group, 'Hybrid')
  }
  else
    group <- c(group, 'Inbred')
}
View(group)

blups <- blups %>% mutate(Group = group, .before= note)
blups <- blups %>% filter(Rep == 1) #removing rep1
View(blups)



## Adding male and female columns for hybrids and inbreds (only female)
male <- c()
female <- c()
for (i in 1:length(blups$genotype)){
  id <- blups$genotype[i]
  if(grepl(' X ', id) == TRUE){
    f <- strsplit(as.character(id), ' X ')[[1]][1]
    m <- strsplit(as.character(id), ' X ')[[1]][2]
    female <- c(female, f)
    male <- c(male, m)
  }
  else {
    male <- c(male, NA)
    female <- c(female, id)
  }
    
}

## Adding SS and NSS information for inbreds


blups <- blups %>% mutate(het_grp = NA, .before=note) # adding a new column for heterotic groups 

hybrids <- blups %>% filter(grepl(' X ', genotype)) #applying a filter for the genotypes containing ' X '

for (i in 1:length(hybrids$genotype)){
  id <- hybrids$genotype[i]
  inbred <- strsplit(as.character(id), ' X ')[[1]][1]
  if(strsplit(as.character(id), ' X ')[[1]][2] == 'B73'){
    blups[blups$genotype == inbred, "het_grp"] <- 'NSS'
  }
  else if(strsplit(as.character(id), ' X ')[[1]][2] == 'Mo17'){
    blups[blups$genotype == inbred, "het_grp"] <- 'SS'
  }
  else if(id == 'B73'){
    blups[blups$genotype == id, "het_grp"] <- 'NSS'
  }
  else {
    blups[blups$genotype == inbred, "het_grp"] <- 'other'
  }
}


for (i in 1:length(blups$genotype)){
  id <- blups$genotype[i]
  if(id == 'B73'){
    blups[blups$genotype == id, "het_grp"] <- 'B73'
  }
  else if (id %in% c('Mo17', 'MO17')){
    blups[blups$genotype == id, "het_grp"] <- 'Mo17'
  }
  else{
    print('empty')
  }
}


res.pca <- prcomp(blups[blups$Group == 'Inbred' & !is.na(blups$het_grp) , 14:length(colnames(blups))], scale = FALSE)

basic_plot <- fviz_pca_ind(res.pca, label= 'none')

ggplot(cbind(basic_plot$data, blups[blups$Group=='Inbred' & !is.na(blups$het_grp), c('Trt', 'Group', 'het_grp')]), aes(x=x, y=y, shape=Trt, col = het_grp ))+
  geom_point(size=2)+
  labs(title = 'PCA', x='Dim1 (82.7%)', y= 'Dim2 (12.4%)')+
  stat_ellipse()+
  theme_bw(14)




blups <- blups %>% mutate(famale = female, male=male, .before = note)

res.pca <- prcomp(blups[blups$Group== 'Hybrid',15:length(colnames(blups))], scale = FALSE)

basic_plot <- fviz_pca_ind(res.pca, label= 'none')

ggplot(cbind(basic_plot$data, blups[blups$Group== 'Hybrid',c('Trt', 'Group', 'male')]), aes(x=x, y=y, shape=Group, col = male ))+
  geom_point(size=2)+
  labs(title = 'PCA', x='Dim1 (82.7%)', y= 'Dim2 (12.4%)')+
  #stat_ellipse(type= 't')+
  theme_bw(14)



View(basic_plot$data)


fviz_pca_ind(res.pca, geom="point")

fviz_pca_ind(res.pca, col.ind="cos2", geom = "point") +
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=0.6)+ theme_minimal()

fviz_pca_ind(res.pca, label="none", habillage=blups[, c('Trt', 'note')])

p <- fviz_pca_ind(res.pca, label="none", habillage=blups[blups$Trt == 'LN', 'note'],
                  addEllipses=TRUE, ellipse.level=0.95)
print(p)

# Control the transparency of variables using their contributions
fviz_pca_var(res.pca, alpha.var="contrib") +
  theme_minimal()


# Select the top 3 contributing variables
fviz_pca_var(res.pca, select.var = list(contrib = 20))


# Change the color by groups, add ellipses
fviz_pca_biplot(res.pca, label= FALSE, habillage=blups[blups$Trt == 'LN', 'note'],
                addEllipses=TRUE, ellipse.level=0.95)


# Use points and text
fviz_pca_var(res.pca, select.var = list(contrib =1000), geom = 'point', pointsize= 0.1)
# Use lines
fviz_pca_var(res.pca, select.var = list(contrib =1000), pointsize= 0.1)


# Select and visualize variables with cos2 >= 0.96
fviz_pca_var(res.pca, select.var = list(cos2 = 0.00006))

              #### PCA of Raw Spectra #####
spectra <- read.csv("Raw_spectrum_merged")

str(spectra)

spectra$Rep <- factor(spectra$Rep)
spectra$Block <- factor(spectra$Block)
spectra$year <- factor(spectra$year)
spectra$genotype <- factor(spectra$genotype)
spectra$note <- factor(spectra$note)
spectra$Trt <- factor(spectra$Trt)
spectra$ASD <- factor(spectra$ASD)

spectra2  <-  spectra %>% drop_na()
spectra <- subset(spectra, select = -c(X, Unnamed..0))

res.pca <- prcomp(spectra2[spectra2$ASD %in% c('1','2'), 13:length(colnames(spectra2))], scale = FALSE)

fviz_pca_ind(res.pca, geom="point")


fviz_pca_ind(res.pca, label="none", habillage=spectra2$ASD)


basic_plot <- fviz_pca_ind(res.pca, label= 'none')
basic_plot
ggplot(cbind(basic_plot$data, spectra2[spectra2$ASD %in% c('1','2'), c('Trt', 'ASD', 'Group')]), aes(x=x, y=y, col=ASD))+
  geom_point(size=1.4, alpha=0.7)+
  labs(title = 'PCA of raw spectrum', x='Dim1 (77.2%)', y= 'Dim2 (10.9%)')+
  #stat_ellipse()+
  theme_bw(14)





                            ###### PCA of blues ######

blues <- read.csv('spectra_blues.csv')
blues <- as.data.frame(blues)
blues <- blues %>% select(-contains(c('.x', '.y', '.z')))
                          
blues <- blues[ , colSums(is.na(blues))==0]


str(blues)


blues$genotype <- as.factor(blues$genotype)
blues$PLOT.ID <- as.factor(blues$PLOT.ID)
blues$rows <- as.factor(blues$rows)
blues$ranges <- as.factor(blues$ranges)
blues$Block <- as.factor(blues$Block)
blues$Rep <- as.factor(blues$Rep)
blues$Trt <- as.factor(blues$Trt)
blues$year <- as.factor(blues$year)
blues$note <-as.factor(blues$note)
blues$Calibration <- as.factor(blues$Calibration)
blues$ASD <- as.factor(blues$ASD)

## adding 'hybrid or inbred' column
group <- c()
for(i in 1:length(blues$note)){
  note <- blues$note[i]
  if(note == 'Hybrid'){
    group <- c(group, 'Hybrid')
  }
  else
    group <- c(group, 'Inbred')
}
View(group)

blues <- blues %>% mutate(Group = group, .before= note)
#blues <- blues %>% filter(Rep == 1) #removing rep1
View(blues)



## Adding male and female columns for hybrids and inbreds (only female)
male <- c()
female <- c()
for (i in 1:length(blues$genotype)){
  id <- blues$genotype[i]
  if(grepl(' X ', id) == TRUE){
    f <- strsplit(as.character(id), ' X ')[[1]][1]
    m <- strsplit(as.character(id), ' X ')[[1]][2]
    female <- c(female, f)
    male <- c(male, m)
  }
  else {
    male <- c(male, NA)
    female <- c(female, id)
  }
  
}

## Adding SS and NSS information for inbreds


blues <- blues %>% mutate(het_grp = NA, .before=note) # adding a new column for heterotic groups 

hybrids <- blues %>% filter(grepl(' X ', genotype)) #applying a filter for the genotypes containing ' X '

for (i in 1:length(hybrids$genotype)){
  id <- hybrids$genotype[i]
  inbred <- strsplit(as.character(id), ' X ')[[1]][1]
  if(strsplit(as.character(id), ' X ')[[1]][2] == 'B73'){
    blues[blues$genotype == inbred, "het_grp"] <- 'NSS'
  }
  else if(strsplit(as.character(id), ' X ')[[1]][2] == 'Mo17'){
    blues[blues$genotype == inbred, "het_grp"] <- 'SS'
  }
  else if(id == 'B73'){
    blues[blues$genotype == id, "het_grp"] <- 'NSS'
  }
  else {
    blues[blues$genotype == inbred, "het_grp"] <- 'other'
  }
}


for (i in 1:length(blues$genotype)){
  id <- blues$genotype[i]
  if(id == 'B73'){
    blues[blues$genotype == id, "het_grp"] <- 'B73'
  }
  else if (id %in% c('Mo17', 'MO17')){
    blues[blues$genotype == id, "het_grp"] <- 'Mo17'
  }
  else{
    print('empty')
  }
}


res.pca <- prcomp(blues[blues$Trt == 'LN' , 11:length(colnames(blues))], scale = FALSE)

basic_plot <- fviz_pca_ind(res.pca, label= 'none')
basic_plot
ggplot(cbind(basic_plot$data, blues[blues$Trt == 'LN', c('Trt', 'note', 'Group')]), aes(x=x, y=y, col=Group))+
  geom_point(size=2)+
  labs(title = 'PCA', x='Dim1 (82.7%)', y= 'Dim2 (12.4%)')+
  #stat_ellipse()+
  theme_bw(14)
 

                                           ##################### PCA of blupsv2  ###########################

blups_merged_v2 <- read.csv('spectra_comb_blups_v2.csv')

str(blups_merged_v2)

blups_merged_v2$genotype <- factor(blups_merged_v2$genotype)
blups_merged_v2$PLOT.ID <- factor(blups_merged_v2$PLOT.ID)
blups_merged_v2$rows <- factor(blups_merged_v2$rows)
blups_merged_v2$ranges <- factor(blups_merged_v2$ranges)
blups_merged_v2$Block <- factor(blups_merged_v2$Block)
blups_merged_v2$Rep <- factor(blups_merged_v2$Rep)
blups_merged_v2$Trt  <- factor(blups_merged_v2$Trt)
blups_merged_v2$year <- factor(blups_merged_v2$year)
blups_merged_v2$note <- factor(blups_merged_v2$note)
blups_merged_v2$Calibration <- factor(blups_merged_v2$Calibration)
blups_merged_v2$ASD <- factor(blups_merged_v2$ASD)
blups_merged_v2$new_GID <- factor(blups_merged_v2$new_GID)

blups_merged_v2 <- blups_merged_v2[blups_merged_v2$Rep == 2 , ]

res.pca <- prcomp(blups_merged_v2[ , 14:length(colnames(blups_merged_v2))], scale = FALSE)

basic_plot <- fviz_pca_ind(res.pca, label= 'none')
basic_plot
ggplot(cbind(basic_plot$data, blups_merged_v2[, c('Trt', 'note', 'Group')]), aes(x=x, y=y, color=Trt))+
  geom_point(aes(shape=Group), size=2)+
  scale_shape_manual(values = c(1,16))+
  labs(title = 'PCA', x='Dim1 (73.8%)', y= 'Dim2 (19.8%)')+
  #stat_ellipse()+
  theme_bw(14)

ind <- get_pca_ind(res.pca)

view(head(ind$contrib))


        ####### Extracting PCs by using Nikee's code ######


color <- blups_merged_v2 ##It should be not have any Nas. This is my main dataframe.

PCA <- prcomp(color[,14:length(colnames(blups_merged_v2))]) ##Passing three columns of red, blue and green values. princomp is inbuilt fucntion.
#color$PC1=PCA$scores[,1] ## Taking first principle component
#color$PC2=PCA$scores[,2] ## Taking second principle component
#color$PC3=PCA$scores[,3] ## Taking third principle component

ind <- get_pca_ind(PCA)

PCs <- ind$coord[, 1:3]

PC.for.gwas <- cbind(color[, c('genotype', "Trt")], PCs)

write.csv(PC.for.gwas, 'GWAS_PCA_blups_v2.csv', row.names = FALSE)

View(PC.for.gwas)

              ###### Extracting PCs from first version of blups ######


color <- read_csv('spectra_blups.csv') ##It should be not have any Nas. This is my main dataframe. 
color <- data_frame(color)
color <- color[which(color$Rep == 1),]
PCA <- prcomp(color[,13:length(colnames(color))]) ## length(colnames(color)) Passing three columns of red, blue and green values. princomp is inbuilt fucntion.
#PCA$scores[,1]
#color$PC1=PCA$scores[,1] ## Taking first principle component
#color$PC2=PCA$scores[,2] ## Taking second principle component
#color$PC3=PCA$scores[,3] ## Taking third principle component

ind <- get_pca_ind(PCA)

PCs <- ind$coord[, 1:3]

PC.for.gwas_blups_v1 <- cbind(color[, c('genotype', "Trt")], PCs)

write.csv(PC.for.gwas_blups_v1, 'PC_GWAS_blups_v1.csv', row.names = FALSE)

View(PC.for.gwas_blups_v1)
