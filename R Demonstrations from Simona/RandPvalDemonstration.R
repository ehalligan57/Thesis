#Simona Dalin
#11/10/20
#Demonstrating linear regression

library(ggplot2)

######## Capture the folder that the current script is in, and change the working directory to this folder so we can save plots easily ############
currentDirName <- dirname(parent.frame(2)$ofile)
print(currentDirName)
setwd(as.character(currentDirName))


########## Arrange The Data ########## 

#This is the data for the plot I sent you
geneExpressionValuesSensitiveCellLine1 <- c(50,43,55,57)
geneExpressionValuesSensitiveCellLine2 <- c(60,53,65,67)
geneExpressionValuesSensitiveCellLine3 <- c(70,63,75,77)
geneExpressionValuesSensitiveCellLine4 <- c(80,73,85,87)

geneExpressionValuesResistantCellLine1 <- c(500,430,550,570)
geneExpressionValuesResistantCellLine2 <- c(600,530,650,670)
geneExpressionValuesResistantCellLine3 <- c(700,630,750,770)
geneExpressionValuesResistantCellLine4 <- c(800,730,850,870)

percentControlGrowthSensitiveCellLine1 <- 10
percentControlGrowthSensitiveCellLine2 <- 17
percentControlGrowthSensitiveCellLine3 <- 20
percentControlGrowthSensitiveCellLine4 <- 30

percentControlGrowthResistantCellLine1 <- 91
percentControlGrowthResistantCellLine2 <- 90
percentControlGrowthResistantCellLine3 <- 89
percentControlGrowthResistantCellLine4 <- 85

#If you run this chunk instead, it will generate random data (and should always give r close to 0, p > 0.05)
geneExpressionValuesSensitiveCellLine1 <- runif(4, min = 0, max = 500)
geneExpressionValuesSensitiveCellLine2 <- runif(4, min = 0, max = 500)
geneExpressionValuesSensitiveCellLine3 <- runif(4, min = 0, max = 500)
geneExpressionValuesSensitiveCellLine4 <- runif(4, min = 0, max = 500)

geneExpressionValuesResistantCellLine1 <- runif(4, min = 0, max = 500)
geneExpressionValuesResistantCellLine2 <- runif(4, min = 0, max = 500)
geneExpressionValuesResistantCellLine3 <- runif(4, min = 0, max = 500)
geneExpressionValuesResistantCellLine4 <- runif(4, min = 0, max = 500)

percentControlGrowthSensitiveCellLine1 <- runif(1, min = 0, max = 100)
percentControlGrowthSensitiveCellLine2 <- runif(1, min = 0, max = 100)
percentControlGrowthSensitiveCellLine3 <- runif(1, min = 0, max = 100)
percentControlGrowthSensitiveCellLine4 <- runif(1, min = 0, max = 100)

percentControlGrowthResistantCellLine1 <- runif(1, min = 0, max = 100)
percentControlGrowthResistantCellLine2 <- runif(1, min = 0, max = 100)
percentControlGrowthResistantCellLine3 <- runif(1, min = 0, max = 100)
percentControlGrowthResistantCellLine4 <- runif(1, min = 0, max = 100)

sensitiveCellLine1 <- data.frame("cellLineName" = "Sensitive1", "geneExpression" = geneExpressionValuesSensitiveCellLine1, "percentControlGrowth" = percentControlGrowthSensitiveCellLine1)
sensitiveCellLine2 <- data.frame("cellLineName" = "Sensitive2", "geneExpression" = geneExpressionValuesSensitiveCellLine2, "percentControlGrowth" = percentControlGrowthSensitiveCellLine2)
sensitiveCellLine3 <- data.frame("cellLineName" = "Sensitive3", "geneExpression" = geneExpressionValuesSensitiveCellLine3, "percentControlGrowth" = percentControlGrowthSensitiveCellLine3)
sensitiveCellLine4 <- data.frame("cellLineName" = "Sensitive4", "geneExpression" = geneExpressionValuesSensitiveCellLine4, "percentControlGrowth" = percentControlGrowthSensitiveCellLine4)

resistantCellLine1 <- data.frame("cellLineName" = "Resistant1", "geneExpression" = geneExpressionValuesResistantCellLine1, "percentControlGrowth" = percentControlGrowthResistantCellLine1)
resistantCellLine2 <- data.frame("cellLineName" = "Resistant2", "geneExpression" = geneExpressionValuesResistantCellLine2, "percentControlGrowth" = percentControlGrowthResistantCellLine2)
resistantCellLine3 <- data.frame("cellLineName" = "Resistant3", "geneExpression" = geneExpressionValuesResistantCellLine3, "percentControlGrowth" = percentControlGrowthResistantCellLine3)
resistantCellLine4 <- data.frame("cellLineName" = "Resistant4", "geneExpression" = geneExpressionValuesResistantCellLine4, "percentControlGrowth" = percentControlGrowthResistantCellLine4)

cellLineData <- rbind(sensitiveCellLine1, sensitiveCellLine2, sensitiveCellLine3, sensitiveCellLine4, resistantCellLine1, resistantCellLine2, resistantCellLine3, resistantCellLine4)




########## Perform the linear regression ########## 
linearRegression <- lm(percentControlGrowth ~ geneExpression, data = cellLineData)
linearRegressionSummary <- summary(linearRegression)
linearRegressionrValue <- sqrt(linearRegressionSummary$r.squared)
linearRegressionfStat <- linearRegressionSummary$fstatistic
linearRegressionPval <- pf(linearRegressionfStat[1],linearRegressionfStat[2],linearRegressionfStat[3],lower.tail=F)

intercept <- linearRegression$coefficients[1]
slope <- linearRegression$coefficients[2]




########## Plot the Data ########## 
plot4 <- ggplot(cellLineData, aes(x = geneExpression, y = percentControlGrowth, color = cellLineName)) + 
  geom_point() +
  geom_abline(intercept = intercept, slope = slope) +
  labs(title = sprintf("r value = %.3s, P-val = %.2e", linearRegressionrValue, linearRegressionPval)) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))


ggsave("Large r and large pval.pdf", plot4, device = "pdf")

