#Caitlin Halligan
#Analysis of Fry et al. Paper
#Genomic predictors of interindividual differences in response to DNA damaging agents 

########## Load Packages ########## 
#import data.table
library(data.table)

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

# Make table with gene expression values of highest and lowest sensitivity groups, treated and untreated (means of two trials)
# Cell lines 6, 4, 9, 20 highest while Cell Lines 12, 8, 22, 7 lowest
TrainingPopulation <- data.frame(Gene=FilteredExpressionValues[,106], CellLine6 = rowMeans(FilteredExpressionValues[,21:22]),CellLine4 = rowMeans(FilteredExpressionValues[,13:14]),CellLine9 = rowMeans(FilteredExpressionValues[,34:35]),CellLine20 = rowMeans(FilteredExpressionValues[,78:79]),CellLine6T=rowMeans(FilteredExpressionValues[,23:24]),CellLine4T = rowMeans(FilteredExpressionValues[,15:16]),CellLine9T = rowMeans(FilteredExpressionValues[,36:37]),CellLine20T = rowMeans(FilteredExpressionValues[,80:81]),CellLine12 = rowMeans(FilteredExpressionValues[,46:47]),CellLine8 = rowMeans(FilteredExpressionValues[,29:30]),CellLine22 = rowMeans(FilteredExpressionValues[,86:87]),CellLine7 = rowMeans(FilteredExpressionValues[,25:26]),CellLine12T = rowMeans(FilteredExpressionValues[,48:49]),CellLine8T = rowMeans(FilteredExpressionValues[,31:32]),CellLine22T = rowMeans(FilteredExpressionValues[,88:89]),CellLine7T = rowMeans(FilteredExpressionValues[,27:28]))

#Find fold difference between training population for basal gene expression
AverageHigh <- rowMeans(TrainingPopulation[,2:5])
AverageLow <- rowMeans(TrainingPopulation[,9:12])
FoldDifference <- data.frame(Gene=FilteredExpressionValues[,106],FoldDifference = (AverageHigh/AverageLow))
#TRUE if within correct range
FoldDifference$new<- (FoldDifference$FoldDifference <= -1.5 & FoldDifference$FoldDifference >= 1.5)
