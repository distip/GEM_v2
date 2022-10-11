library(MASS)
library(ggplot2)


blups <- read.csv('spectra_blups.csv')
blups<- blups[blups$Rep == 1, ]

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
blups$Group <- as.factor(blups$Group)


blups.sub <- subset(blups, select = -c(genotype, PLOT.ID, rows, ranges, Block, Rep, Trt, year,Group, Calibration, ASD))

#Use 70% of dataset as training set and remaining 30% as testing set
sample <- sample(c(TRUE, FALSE), nrow(blups.sub), replace=TRUE, prob=c(0.7,0.3))
train <- blups.sub[sample, ]
test <- blups.sub[!sample, ] 

#fit LDA model
model <- lda(note~., data=train)

#view model output
model[1]


#use LDA model to make predictions on test data
predicted <- predict(model, test)

names(predicted)


#view predicted class for first six observations in test set
head(predicted$class)


#view posterior probabilities for first six observations in test set
head(predicted$posterior)


#view linear discriminants for first six observations in test set
head(predicted$x)


#find accuracy of model
mean(predicted$class==test$note)

#define data to plot
lda_plot <- cbind(train, predict(model)$x)

#create plot
ggplot(lda_plot, aes(LD1, LD2)) +
  geom_point(aes(color = note), size=5)
