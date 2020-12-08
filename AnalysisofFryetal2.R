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

#This line will set the working directory to the folder that this script is saved within
setwd(dirname(sys.frame(1)$ofile))

#/Users/caitlinhalligan (Home folder) is working directory, call that dirName
dirName <- getwd() 

########## Load Data ########## 
#Load data
#ProbsetExpressionData <- fread("/Users/caitlinhalligan/GSE10313_series_matrix.txt",skip="ID_REF",nrows=54613)
ProbsetExpressionData <- fread("./data/GSE10313_series_matrix.txt",skip="ID_REF",nrows=54613)


########## Add Gene Names to the Data ########## 

#Load GeneNames Data
#ProbsettoGene <- fread("/Users/caitlinhalligan/GeneNames.csv",select=c("Probe Set ID","Gene Symbol"),skip="Probe Set ID")
ProbsettoGene <- fread("./data/GeneNames.csv",select=c("Probe Set ID","Gene Symbol"),skip="Probe Set ID")


#If the values of ID_REF and Probe Set ID are equal, then change ID_REF to Gene Symbol
ProbsetExpressionData$GeneSymbol <- NA
for(i in 1:nrow(FilteredExpressionValues)){
  rowsCorrespondingToThisProbeSet <- ProbsetExpressionData$ID_REF[i] == ProbsettoGene$`Probe Set ID`
  ProbsetExpressionData$GeneSymbol[i] <- ProbsettoGene$`Gene Symbol`[rowsCorrespondingToThisProbeSet]
}  


#Same as that for loop, but using sapply
findGeneSymbol <- function(probesetID){
  rowsCorrespondingToThisProbeSet <- probesetID == ProbsettoGene$`Probe Set ID`
  ProbsettoGene$`Gene Symbol`[rowsCorrespondingToThisProbeSet]
}
geneSymbols <- sapply(ProbsetExpressionData$ID_REF, findGeneSymbol)
ProbsetExpressionData$GeneSymbol <- geneSymbols


########## Remove all rows corresponding to control probes/etc (ie, remove all rows with a blank gene symbol) ############
blankRows <- ProbsetExpressionData$GeneSymbol == "" 
gapRows <- ProbsetExpressionData$GeneSymbol == "---" 
rowsToRemove <- blankRows | gapRows
ProbsetExpressionDataNonGenesRemoved <- ProbsetExpressionData[!rowsToRemove,] 


########## Filter out Low Expressed Probesets ########## 
#The paper says we should end up with 19290 probesets, so we'll set a cutoff
#At the mean intensity across arrays that gives 19290 probesets above that value

#Mean of every row and every column except 1st (which has probeset names)
#rowMeans(ProbsetExpressionData[,2:ncol(ProbsetExpressionData)])
#excludes 1st column
MeanofRows<-rowMeans(ProbsetExpressionDataNonGenesRemoved[,2:105])
thresholdBelow <- sort(MeanofRows,decreasing=TRUE)[19290]
thresholdAbove <- sort(MeanofRows,decreasing=TRUE)[19291]
threshold<-mean(c(thresholdAbove,thresholdBelow))
#only have rows above threshold 
FilteredExpressionValues<-ProbsetExpressionDataNonGenesRemoved[MeanofRows>threshold,]




########## Data Table containing Training Population ###########

# Make table with gene expression values of highest and lowest sensitivity groups, treated and untreated (means of two trials)
# Cell lines 6, 4, 9, 20 highest while Cell Lines 12, 8, 22, 7 lowest

#This information is from the GEO website for this accession number: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10313
CellLine6Columns <- colnames(FilteredExpressionValues) %in% c("GSM260573","GSM260574")
CellLine4Columns <- colnames(FilteredExpressionValues) %in% c("GSM260565","GSM260566")
CellLine9Columns <- colnames(FilteredExpressionValues) %in% c("GSM260585","GSM260586")
CellLine20Columns <- colnames(FilteredExpressionValues) %in% c("GSM260629","GSM260630")

CellLine6TColumns <- colnames(FilteredExpressionValues) %in% c("GSM260575","GSM260576")
CellLine4TColumns <- colnames(FilteredExpressionValues) %in% c("GSM260567","GSM260568")
CellLine9TColumns <- colnames(FilteredExpressionValues) %in% c("GSM260587","GSM260588")
CellLine20TColumns <- colnames(FilteredExpressionValues) %in% c("GSM260631","GSM260632")

CellLine12Columns <- colnames(FilteredExpressionValues) %in% c("GSM260597","GSM260598")
CellLine8Columns <- colnames(FilteredExpressionValues) %in% c("GSM260581","GSM260582")
CellLine22Columns <- colnames(FilteredExpressionValues) %in% c("GSM260637","GSM260638")
CellLine7Columns <- colnames(FilteredExpressionValues) %in% c("GSM260577","GSM260578")

CellLine12TColumns <- colnames(FilteredExpressionValues) %in% c("GSM260599","GSM260600")
CellLine8TColumns <- colnames(FilteredExpressionValues) %in% c("GSM260583","GSM260584")
CellLine22TColumns <- colnames(FilteredExpressionValues) %in% c("GSM260639","GSM260574")
CellLine7TColumns <- colnames(FilteredExpressionValues) %in% c("GSM260579","GSM260580")

#Use those columns to extract the training dataset
TrainingPopulation <- data.frame(Gene = FilteredExpressionValues$GeneSymbol, 
                                 probeSet = FilteredExpressionValues$ID_REF, 
                                 CellLine6 = rowMeans(FilteredExpressionValues[, ..CellLine6Columns]),
                                 CellLine4 = rowMeans(FilteredExpressionValues[, ..CellLine4Columns]),
                                 CellLine9 = rowMeans(FilteredExpressionValues[, ..CellLine9Columns]),
                                 CellLine20 = rowMeans(FilteredExpressionValues[, ..CellLine20Columns]),
                                 CellLine6T = rowMeans(FilteredExpressionValues[, ..CellLine6TColumns]),
                                 CellLine4T = rowMeans(FilteredExpressionValues[, ..CellLine4TColumns]),
                                 CellLine9T = rowMeans(FilteredExpressionValues[, ..CellLine9TColumns]),
                                 CellLine20T = rowMeans(FilteredExpressionValues[, ..CellLine20TColumns]),
                                 CellLine12 = rowMeans(FilteredExpressionValues[, ..CellLine12Columns]),
                                 CellLine8 = rowMeans(FilteredExpressionValues[, ..CellLine8Columns]),
                                 CellLine22 = rowMeans(FilteredExpressionValues[, ..CellLine22Columns]),
                                 CellLine7 = rowMeans(FilteredExpressionValues[, ..CellLine7Columns]),
                                 CellLine12T = rowMeans(FilteredExpressionValues[, ..CellLine12TColumns]),
                                 CellLine8T = rowMeans(FilteredExpressionValues[, ..CellLine8TColumns]),
                                 CellLine22T = rowMeans(FilteredExpressionValues[, ..CellLine22TColumns]),
                                 CellLine7T = rowMeans(FilteredExpressionValues[, ..CellLine7TColumns]))

########## Fold Difference, TtestPvalue, LinearRegressionR and Pval for Unique Genes ##########
UniqueGeneSymbol <- unique(FilteredExpressionValues$GeneSymbol)
TrainingDataFilteringValues<-data.frame(Gene=TrainingPopulation$Gene,
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

#Set which columns of TrainingPopulation correspond to resistant vs. sensitive lines


# for (GeneIndex in 1:length(UniqueGeneSymbol)){
#   LogicalUniqueGenes<- TrainingPopulation$Gene==UniqueGeneSymbol[GeneIndex]
#   AverageHigh <- rowMeans(as.matrix(TrainingPopulation[LogicalUniqueGenes,2:5]))
#   AverageLow <- rowMeans(as.matrix(TrainingPopulation[LogicalUniqueGenes,10:13]))
#   TrainingDataFilteringValues$FoldDifference[GeneIndex] <- log2(AverageHigh/AverageLow)
#   TrainingDataFilteringValues$TtestPvalue[GeneIndex] <- (t.test(TrainingPopulation[LogicalUniqueGenes,2:5],TrainingPopulation[LogicalUniqueGenes,10:13], var.equal = FALSE))$p.value
#   
#   numberOfProbesetsCorrespondingToThisGene <- dim(TrainingPopulation[LogicalUniqueGenes,c(2:5,10:13)])[1]
#   sensitivityAndExpressionDataThisGene <- data.frame("percentControlGrowth" = rep(cellLineData$percentControlGrowth,numberOfProbesetsCorrespondingToThisGene),
#                                                      "geneExpression" = array(t(TrainingPopulation[LogicalUniqueGenes,c(2:5,10:13)])))
#   linearRegression <- lm(geneExpression ~ percentControlGrowth, sensitivityAndExpressionDataThisGene)
#   linearRegressionSummary <- summary(linearRegression)
#   linearRegressionfStat <- linearRegressionSummary$fstatistic
#   TrainingDataFilteringValues$LinearRegressionR[GeneIndex] <- sqrt(linearRegressionSummary$r.squared)
#   TrainingDataFilteringValues$LinearRegressionPval[GeneIndex] <- pf(linearRegressionfStat[1],linearRegressionfStat[2],linearRegressionfStat[3],lower.tail=F)
# }

#Calculating fold change and correlation by probe rather than by gene
for (probesetIndex in 1:(dim(TrainingPopulation)[1])){
  AverageHigh <- mean(as.matrix(TrainingPopulation[probesetIndex,3:6]))
  AverageLow <- mean(as.matrix(TrainingPopulation[probesetIndex,11:14]))
  TrainingDataFilteringValues$FoldDifference[probesetIndex] <- log2(AverageHigh/AverageLow)
  TrainingDataFilteringValues$TtestPvalue[probesetIndex] <- (t.test(TrainingPopulation[probesetIndex,3:6],TrainingPopulation[probesetIndex,11:14], var.equal = FALSE))$p.value
  
  numberOfProbesetsCorrespondingToThisGene <- dim(TrainingPopulation[probesetIndex,c(3:6,11:14)])[1]
  sensitivityAndExpressionDataThisGene <- data.frame("percentControlGrowth" = rep(cellLineData$percentControlGrowth,numberOfProbesetsCorrespondingToThisGene),
                                                     "geneExpression" = array(t(TrainingPopulation[probesetIndex,c(3:6,11:14)])))
  correlation <- cor.test(sensitivityAndExpressionDataThisGene$percentControlGrowth,sensitivityAndExpressionDataThisGene$geneExpression,method = "pearson")
  TrainingDataFilteringValues$LinearRegressionR[probesetIndex] <- correlation$estimate
  TrainingDataFilteringValues$LinearRegressionPval[probesetIndex] <- correlation$p.value
}

#Add column which is TRUE if gene is in BASA set
TrainingDataFilteringValues$BASA <- abs(TrainingDataFilteringValues$FoldDifference) >= 1.5 & 
  TrainingDataFilteringValues$TtestPvalue < 0.05 & 
  abs(TrainingDataFilteringValues$LinearRegressionR) >= 0.7 & 
  TrainingDataFilteringValues$LinearRegressionPval < 0.01

TrainingDataFilteringValues$differentialExpression <- abs(TrainingDataFilteringValues$FoldDifference) >= 1.5 & 
  TrainingDataFilteringValues$TtestPvalue < 0.05

TrainingDataFilteringValues$correlation <- abs(TrainingDataFilteringValues$LinearRegressionR) >= 0.7 & 
  TrainingDataFilteringValues$LinearRegressionPval < 0.01



########## Plot the Data ########## 
cellLineData2 <- data.frame("percentControlGrowth" = cellLineData$percentControlGrowth, "geneExpression" = as.numeric(TrainingPopulation[TrainingPopulation$Gene == "MGMT",2:13][,-5:-8]), "cellLineName" = colnames(TrainingPopulation)[2:13][-5:-8])

LinearRegressionGraph <- ggplot(cellLineData2, aes(x = percentControlGrowth, y = geneExpression, color = cellLineName)) + 
  geom_point() +
  geom_abline(intercept = intercept, slope = slope) +
  labs(title = sprintf("r value = %.3s, P-val = %.2e", linearRegressionrValue, linearRegressionPval)) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))


ggsave("LinearRegressionGraph.pdf", LinearRegressionGraph, device = "pdf")


