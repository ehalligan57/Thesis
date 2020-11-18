#Caitlin Halligan
#Analysis of Fry et al. Paper
#Genomic predictors of interindividual differences in response to DNA damaging agents 

########## Load Packages ########## 
#import data.table
library(data.table)
library(tidyverse)

#Install ggplot2
#install.packages("ggplot2")
library(ggplot2)

#/Users/caitlinhalligan (Home folder) is working directory, call that dirName
dirName <- getwd() 

########## Load Data ########## 
#Load data
ProbsetExpressionData <- fread("/Users/caitlinhalligan/GSE10313_series_matrix.txt",skip="ID_REF",nrows=54613)

########## Filter out Low Expressed Probesets ########## 
#The paper says we should end up with 19290 probesets, so we'll set a cutoff
#At the mean intensity across arrays that gives 19290 probesets above that value

#Mean of every row and every column except 1st (which has probeset names)
rowMeans(ProbsetExpressionData[,2:ncol(ProbsetExpressionData)])
#excludes 1st column
MeanofRows<-rowMeans(ProbsetExpressionData[,-1])
thresholdBelow <- sort(MeanofRows,decreasing=TRUE)[19290]
thresholdAbove <- sort(MeanofRows,decreasing=TRUE)[19291]
threshold<-mean(c(thresholdAbove,thresholdBelow))
#only have rows above threshold 
FilteredExpressionValues<-ProbsetExpressionData[MeanofRows>threshold,]


########## Add Gene Names to the Filtered Data ########## 

#Load GeneNames Data
ProbsettoGene <- fread("/Users/caitlinhalligan/GeneNames.csv",select=c("Probe Set ID","Gene Symbol"),skip="Probe Set ID")

#If the values of ID_REF and Probe Set ID are equal, then change ID_REF to Gene Symbol
FilteredExpressionValues$GeneSymbol <- NA
for(i in 1:nrow(FilteredExpressionValues)){
  rowsCorrespondingToThisProbeSet <- FilteredExpressionValues$ID_REF[i] == ProbsettoGene$`Probe Set ID`
  FilteredExpressionValues$GeneSymbol[i] <- ProbsettoGene$`Gene Symbol`[rowsCorrespondingToThisProbeSet]
}

########## Data Table containing Training Population ###########

# Make table with gene expression values of highest and lowest sensitivity groups, treated and untreated (means of two trials)
# Cell lines 6, 4, 9, 20 highest while Cell Lines 12, 8, 22, 7 lowest
TrainingPopulation <- data.frame(Gene=FilteredExpressionValues[,106], 
                                 CellLine6 = rowMeans(FilteredExpressionValues[,21:22]),
                                 CellLine4 = rowMeans(FilteredExpressionValues[,13:14]),
                                 CellLine9 = rowMeans(FilteredExpressionValues[,34:35]),
                                 CellLine20 = rowMeans(FilteredExpressionValues[,78:79]),
                                 CellLine6T=rowMeans(FilteredExpressionValues[,23:24]),
                                 CellLine4T = rowMeans(FilteredExpressionValues[,15:16]),
                                 CellLine9T = rowMeans(FilteredExpressionValues[,36:37]),
                                 CellLine20T = rowMeans(FilteredExpressionValues[,80:81]),
                                 CellLine12 = rowMeans(FilteredExpressionValues[,46:47]),
                                 CellLine8 = rowMeans(FilteredExpressionValues[,29:30]),
                                 CellLine22 = rowMeans(FilteredExpressionValues[,86:87]),
                                 CellLine7 = rowMeans(FilteredExpressionValues[,25:26]),
                                 CellLine12T = rowMeans(FilteredExpressionValues[,48:49]),
                                 CellLine8T = rowMeans(FilteredExpressionValues[,31:32]),
                                 CellLine22T = rowMeans(FilteredExpressionValues[,88:89]),
                                 CellLine7T = rowMeans(FilteredExpressionValues[,27:28]))

########## Fold Difference, TtestPvalue, LinearRegressionR and Pval for Unique Genes ##########
UniqueGeneSymbol <- unique(FilteredExpressionValues$GeneSymbol)
TrainingDataFilteringValues<-data.frame(Gene=UniqueGeneSymbol,
                                        FoldDifference=NaN,
                                        TtestPvalue=NaN,
                                        LinearRegressionR=NaN,
                                        LinearRegressionPval=NaN)
#estimate percent control growth from Fry et al. paper
percentControlGrowthSensitiveCellLine6 <- 19
percentControlGrowthSensitiveCellLine4 <- 22
percentControlGrowthSensitiveCellLine9 <- 35
percentControlGrowthSensitiveCellLine20 <- 36

percentControlGrowthResistantCellLine12 <- 80
percentControlGrowthResistantCellLine8 <- 81
percentControlGrowthResistantCellLine22 <- 83
percentControlGrowthResistantCellLine7 <- 85

sensitiveCellLine6 <- data.frame("cellLineName" = "Sensitive6", "percentControlGrowth" = percentControlGrowthSensitiveCellLine6)
sensitiveCellLine4 <- data.frame("cellLineName" = "Sensitive4", "percentControlGrowth" = percentControlGrowthSensitiveCellLine4)
sensitiveCellLine9 <- data.frame("cellLineName" = "Sensitive9", "percentControlGrowth" = percentControlGrowthSensitiveCellLine9)
sensitiveCellLine20 <- data.frame("cellLineName" = "Sensitive20", "percentControlGrowth" = percentControlGrowthSensitiveCellLine20)

resistantCellLine12 <- data.frame("cellLineName" = "Resistant12", "percentControlGrowth" = percentControlGrowthResistantCellLine12)
resistantCellLine8 <- data.frame("cellLineName" = "Resistant8", "percentControlGrowth" = percentControlGrowthResistantCellLine8)
resistantCellLine22 <- data.frame("cellLineName" = "Resistant22", "percentControlGrowth" = percentControlGrowthResistantCellLine22)
resistantCellLine7 <- data.frame("cellLineName" = "Resistant7", "percentControlGrowth" = percentControlGrowthResistantCellLine7)

cellLineData <- rbind(sensitiveCellLine6, sensitiveCellLine4, sensitiveCellLine9, sensitiveCellLine20, 
                      resistantCellLine12, resistantCellLine8, resistantCellLine22, resistantCellLine7)

for (GeneIndex in 1:length(UniqueGeneSymbol)){
  LogicalUniqueGenes<- TrainingPopulation$Gene==UniqueGeneSymbol[GeneIndex]
  AverageHigh <- mean(as.matrix(TrainingPopulation[LogicalUniqueGenes,2:5]))
  AverageLow <- mean(as.matrix(TrainingPopulation[LogicalUniqueGenes,10:13]))
  TrainingDataFilteringValues$FoldDifference[GeneIndex] <- abs(AverageHigh/AverageLow)
  TrainingDataFilteringValues$TtestPvalue[GeneIndex] <- (t.test(TrainingPopulation[LogicalUniqueGenes,2:5],TrainingPopulation[LogicalUniqueGenes,10:13], var.equal = TRUE))$p.value
  linearRegression <- lm(cellLineData$percentControlGrowth ~ as.numeric(TrainingPopulation[GeneIndex,c(2:5,10:13)]), data = TrainingPopulation)
  linearRegressionSummary <- summary(linearRegression)
  linearRegressionfStat <- linearRegressionSummary$fstatistic
  TrainingDataFilteringValues$LinearRegressionR[GeneIndex] <- abs(sqrt(linearRegressionSummary$r.squared))
  TrainingDataFilteringValues$LinearRegressionPval[GeneIndex] <- pf(linearRegressionfStat[1],linearRegressionfStat[2],linearRegressionfStat[3],lower.tail=F)
}

#Add column which is TRUE if gene is in BASA set
TrainingDataFilteringValues$BASA<-TrainingDataFilteringValues['FoldDifference'] >= 1.5 & TrainingDataFilteringValues['TtestPvalue'] < 0.05 & TrainingDataFilteringValues['LinearRegressionR'] >= 0.7 & TrainingDataFilteringValues['LinearRegressionPval'] < 0.01


########## Plot the Data ########## 
LinearRegressionGraph <- ggplot(cellLineData, aes(x = cellLineData$percentControlGrowth, y = as.numeric(TrainingPopulation[GeneIndex,c(2:5,10:13)]), color = cellLineName)) + 
  geom_point() +
  geom_abline(intercept = intercept, slope = slope) +
  labs(title = sprintf("r value = %.3s, P-val = %.2e", linearRegressionrValue, linearRegressionPval)) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))


ggsave("LinearRegressionGraph.pdf", LinearRegressionGraph, device = "pdf")


