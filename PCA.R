library(FactoMineR)
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

## Adding a new column named 'Group' indicating hybrids and inbreds

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


res.pca <- prcomp(blups[,13:length(colnames(blups))], scale = FALSE)

basic_plot <- fviz_pca_ind(res.pca, label= 'none')

ggplot(cbind(basic_plot$data, blups[,c('Trt', 'Group')]), aes(x=x, y=y, col = Group, shape=Trt))+
  geom_point(size=1.7)+
  labs(title = 'PCA', x='Dim1 (82.7%)', y= 'Dim2 (12.4%)')


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

spectra2  <-  spectra %>% drop_na()

res.pca <- prcomp(spectra2[spectra2$Trt == 'LN', 12:length(colnames(spectra2))], scale = FALSE)

fviz_pca_ind(res.pca, geom="point")


fviz_pca_ind(res.pca, label="none", habillage=spectra2[spectra2$Trt == 'LN', c('ASD')])
