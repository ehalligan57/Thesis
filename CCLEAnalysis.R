#CCLE analysis

########## Load Packages ########## 
library(data.table)
library(tidyverse)
library(ggplot2)

########## Load Data ##########

#Setting working directory for either of our computers. Uncomment the line for your computer
setwd("/Users/sdalin/Dropbox (Partners HealthCare)/Postdoc/Projects/CCNU/Caitlin/Thesis")
#setwd("/Users/caitlinhalligan/Desktop/Thesis")

RNAseqcounts <- fread("./data/CCLE_RNAseq_genes_counts_20180929.gct") 
RNAseqcountsGlioma <- select(RNAseqcounts,Name, Description,contains("CENTRAL_NERVOUS_SYSTEM"))
Methylation <- fread("./data/CCLE_RRBS_TSS1kb_20181022.txt") 
# Contains experiment_id (convert to cell line name) and master_cpd_id (convert to drug name)#
Sensitivity <- fread("./data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.data.curves_post_qc.txt") 
# File has master_ccl_id and ccl_name #
CellLineMeta <- fread("./data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_cell_line.txt") 
# File has master_cpd_id to cdp_name (the drug name) #
CompoundMeta <- fread("./data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_compound.txt") 
# File has experiment_id and master_ccl_id #
ExperimentMeta <- fread("./data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_experiment.txt") 

########## Filter for MGMT Data ##########
MGMTexpression <- as.data.frame(t(RNAseqcounts[which(RNAseqcounts$Description=="MGMT"),])) 
MGMTmethylation <- as.data.frame(t(Methylation[grep("MGMT", Methylation$locus_id),])) 

# Add row to Sensitivity data which has cell line names based on experiment ID
# First, convert experiment_id to master_ccl_id (from ExperimentMeta)
# Next, convert master_ccl_id to cpd_name (CellLineMeta)
Sensitivity$CellLineName <- NA
for(ExperimentID in unique(Sensitivity$experiment_id)){ #loops over each row
  if (length(which(ExperimentID==ExperimentMeta$experiment_id))==0){ 
    cclID <- NA
  }else{
    cclID <- ExperimentMeta$master_ccl_id[ExperimentID==ExperimentMeta$experiment_id]
    ExperimentIDrows <- Sensitivity$experiment_id==ExperimentID
    Sensitivity$CellLineName[ExperimentIDrows] <- CellLineMeta$ccl_name[cclID==CellLineMeta$master_ccl_id]
  } 
}

########## Filter out just alkylators ##########
DrugNametoID <-filter(CompoundMeta,grepl("DNA alkylator",target_or_activity_of_compound)) #cpd_name to master_cpd_id
AlkylatorSensitivity<-filter(Sensitivity, master_cpd_id %in% DrugNametoID$master_cpd_id) 
TMZSensitivity <- filter(Sensitivity, master_cpd_id=="351043") 
#check if it worked:
DrugNametoID$master_cpd_id
unique(AlkylatorSensitivity$master_cpd_id) # these two are the same

plot(Sensitivity$conc_pts_fit,Sensitivity$area_under_curve) 
#there is not a trend in AUC and points fit, so there may be a normalization going on that we don't know about

########## Create Data Frame with MGMT methylation, MGMT expression, TMZ sensitivity for all cell lines ##########
EdgeRdataframe <- data.frame(celllines=colnames(RNAseqcounts)[-(1:2)]) #create table with first column as cell line names
EdgeRdataframe$MGMTexpressiondata <- NA 
for(RowNumber in (1:nrow(EdgeRdataframe))){ #loops over cell line names
  CellLineName <- EdgeRdataframe$celllines[RowNumber] 
  CellLineNumber <- unlist(strsplit(as.character(CellLineName),"_"))[1]
  EdgeRdataframe$MGMTexpressiondata[RowNumber] <- MGMTexpression$V1[CellLineName==row.names(MGMTexpression)]  
}
#for cell lines which have multiple TMZ sensitivity trials, the cell lines could have used different media. A possible
# improvement is using the media which is most commonly used for those cell lines

########## Add Drug Sensitivity to Table ##########
for (AlkylatorcpdID in unique(AlkylatorSensitivity$master_cpd_id)){ #for each Alkylating agent
  Colnumber<-ncol(EdgeRdataframe)
  CurrentAlkylatorSensitivityData <- subset(AlkylatorSensitivity, master_cpd_id == AlkylatorcpdID)
  AlkylatorSensitivityAverage<- aggregate(CurrentAlkylatorSensitivityData$area_under_curve,by=list(CellLineName=CurrentAlkylatorSensitivityData$CellLineName),data=CurrentAlkylatorSensitivityData,FUN=mean)
  for(RowNumber in (1:nrow(EdgeRdataframe))){ #loops over cell line names
    CellLineName <- EdgeRdataframe$celllines[RowNumber] 
    CellLineNumber <- unlist(strsplit(as.character(CellLineName),"_"))[1]
    if (length(which(CellLineNumber==AlkylatorSensitivityAverage$CellLineName))==0){
      EdgeRdataframe[RowNumber,Colnumber+1] <- NA
    }else{
      EdgeRdataframe[RowNumber,Colnumber+1] <- AlkylatorSensitivityAverage$x[which(CellLineNumber==AlkylatorSensitivityAverage$CellLineName)]
    } 
  }
  colnames(EdgeRdataframe)[Colnumber+1]<-paste(DrugNametoID$cpd_name[AlkylatorcpdID==DrugNametoID$master_cpd_id],'SensitivityAUC')
}

########## Filter table for just glioma cell lines ##########
count(grep("CENTRAL_NERVOUS_SYSTEM",EdgeRdataframe$celllines,fixed=TRUE)) #65 glioma cell lines
EdgeRdataframe <- EdgeRdataframe[grep("CENTRAL_NERVOUS_SYSTEM",EdgeRdataframe$celllines,fixed=TRUE),]
# how many useable cell lines for EdgeR? 
65 - sum(is.na(EdgeRdataframe$`temozolomide SensitivityAUC`)) #45

########## Label cell lines as high or low MGMT ##########
ggplot(data = EdgeRdataframe, aes(x = EdgeRdataframe$MGMTexpressiondata)) + geom_histogram() + geom_histogram(binwidth=10)
median(EdgeRdataframe$MGMTexpressiondata) #311
for(RowNumber in (1:nrow(EdgeRdataframe))){ #loops over rows
  if (EdgeRdataframe$MGMTexpressiondata[RowNumber] < 311){
    EdgeRdataframe$MGMTexpressiondata[RowNumber] = 'low'
  }else{
    EdgeRdataframe$MGMTexpressiondata[RowNumber] = 'high'
  }
}
length(which(EdgeRdataframe$MGMTexpressiondata=='high')) #42 with high, 23 with low

########## Label cell lines as High/low sensitivity for alkylating agents ##########
ggplot(data = EdgeRdataframe, aes(x = EdgeRdataframe$`chlorambucil SensitivityAUC`)) + geom_histogram() + geom_histogram(binwidth=.05)
for(ColNumber in (3:ncol(EdgeRdataframe))){
  MedianDrug = median(EdgeRdataframe[,ColNumber],na.rm=TRUE)
  for(RowNumber in (1:nrow(EdgeRdataframe))){
    if (is.na(EdgeRdataframe[RowNumber,ColNumber])){
      EdgeRdataframe[RowNumber,ColNumber] = NA
    }else if (EdgeRdataframe[RowNumber,ColNumber] < MedianDrug){
      EdgeRdataframe[RowNumber,ColNumber] = 'sensitive'
    }else{
      EdgeRdataframe[RowNumber,ColNumber] = 'resistant'
    }
  }
}

########     Plots describing the cell lines     #########
# MGMT expression and methylation
ggplot(data = EdgeRdataframe, aes(x = MGMTmethylationdata, y = MGMTexpressiondata)) + geom_point()
#why is it not linear? #use gene expression and not methylation 

# MGMT methylation or expression vs TMZ
ggplot(data = EdgeRdataframe, aes(x = MGMTmethylationdata, y = `temozolomide SensitivityAUC`)) + geom_point()
ggplot(data = EdgeRdataframe, aes(x = MGMTexpressiondata, y = `temozolomide SensitivityAUC`)) + geom_point()

#MGMT methylation of cell lines
ggplot(data = EdgeRdataframe, aes(x = MGMTmethylationdata)) + geom_histogram() + geom_histogram(binwidth=10)
#many have low MGMT methylation
#Graph theme
theme<-theme(panel.background = element_blank(),
             panel.border=element_rect(color = "black", size = 1, fill=NA),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             axis.text.x=element_text(colour="black", size = 15),
             axis.text.y=element_text(colour="black", size = 15),
             axis.title.y = element_text(colour = "black", size = 15),
             axis.title.x = element_text(colour = "black", size = 15),
             axis.ticks=element_line(colour="black"),
             plot.margin=unit(c(1,1,1,1),"line"),
             legend.title = element_blank(),
             legend.key = element_blank(),
             plot.title = element_text(face = "bold", hjust = 0.5),
             legend.background = element_rect(fill=alpha('blue', 0)),
             legend.text = element_text(size=14))

MGMTExpressionplot <- ggplot(data = EdgeRdataframe, aes(x = MGMTexpressiondata)) + geom_histogram() + geom_histogram(binwidth=1) + theme
print(MGMTExpressionplot + 
        ggtitle("Distribution of MGMT Expression Across Glioma Cell Lines")+
        labs(y="Number of Cell Lines", x = "MGMT Expression (counts)"))

#TMZ sensitivity for the cell lines 
TMZSensitivityplot <- ggplot(data = EdgeRdataframe, aes(x = `temozolomide SensitivityAUC`)) + geom_histogram() + geom_histogram(binwidth=.2)
print(TMZSensitivityplot + 
        ggtitle("Distribution of TMZ Sensitivity Across Glioma Cell Lines")+
        labs(y="Number of Cell Lines", x = "TMZ Sensitivity (AUC values)")) + theme


#########     Making the Design Matrix    ######### 
#make new dataframe with just MGMT high cell lines and just MGMT low cell lines
#for each drug(col) if MGMT in that row is high create designmat separating res and sen and if MGMT is low create designmat separating res and sen
ListofDesignMatricesMGMThigh <- c()
ListofDesignMatricesMGMTlow <- c()
CellLinesMGMTHigh <- c() #Need to initialize these guys
CellLinesMGMTLow <- c() #Need to initialize these guys
for(ColNumber in (3:ncol(EdgeRdataframe))){#loops over drugs
  MGMThighSensitivity <- c()
  MGMTlowSensitivity <- c()
  CellLineNamesMGMThigh <- c() #Need to initialize with the same name as you use in the loop
  CellLineNamesMGMTlow <- c() #Need one for MGMTlow as well
  for(RowNumber in (1:nrow(EdgeRdataframe))){ #loops of cell lines
    if (EdgeRdataframe$MGMTexpressiondata[RowNumber]=='high'){
      if (is.na(EdgeRdataframe[RowNumber,ColNumber])){ #skip NAs
        next
      }else{
        new_sensitivity <- EdgeRdataframe[RowNumber,ColNumber]
        MGMThighSensitivity <- c(MGMThighSensitivity, new_sensitivity)
        new_cellline <- EdgeRdataframe[RowNumber,1] #each run of the for loop it stores the cell line name
        CellLineNamesMGMThigh <- c(CellLineNamesMGMThigh, new_cellline) #all cell line names in vector for MGMT high
      }
    }else{ #MGMT low
      if (is.na(EdgeRdataframe[RowNumber,ColNumber])){ #skip NAs
      }else{
        new_sensitivity <- EdgeRdataframe[RowNumber,ColNumber]
        MGMTlowSensitivity <- c(MGMTlowSensitivity, new_sensitivity)
        new_cellline <- EdgeRdataframe[RowNumber,1] #each run of the for loop it stores the cell line name
        CellLineNamesMGMTlow <- c(CellLineNamesMGMTlow, new_cellline) #all cell line names in vector for MGMT low
      }
    }
  }
  new_designMatMGMThigh <- model.matrix(~0+MGMThighSensitivity) 
  new_designMatMGMTlow <- model.matrix(~0+MGMTlowSensitivity)

  nextEmptySlot <- length(ListofDesignMatricesMGMThigh) + 1
  
  ListofDesignMatricesMGMThigh[[nextEmptySlot]] <- new_designMatMGMThigh
  CellLinesMGMTHigh[[nextEmptySlot]] <- CellLineNamesMGMThigh #list of cell line names for each drug in order of ListofDesignMatricesMGMThigh
  ListofDesignMatricesMGMTlow[[nextEmptySlot]] <- new_designMatMGMTlow
  CellLinesMGMTLow[[nextEmptySlot]] <- CellLineNamesMGMTlow #list of cell line names for each drug in order of ListofDesignMatricesMGMTlow
} 

#Name the slots of the lists of design matrices by the drug they refer to
names(ListofDesignMatricesMGMThigh) <- as.vector(sapply(names(EdgeRdataframe)[3:ncol(EdgeRdataframe)], function(x){substr(x, 1, regexpr(" ", x)-1)}))
names(ListofDesignMatricesMGMTlow) <- sapply(names(EdgeRdataframe)[3:ncol(EdgeRdataframe)], function(x){substr(x, 1, regexpr(" ", x)-1)})


#########     EdgeR Workflow    ######### 
ListofAllDesignMatrices <- c(ListofDesignMatricesMGMThigh,ListofDesignMatricesMGMTlow)
ListofAllCellLines <- c(CellLinesMGMTHigh,CellLinesMGMTLow)
for(DesignMatrixIndex in 1:length(ListofAllDesignMatrices)){
  #trying to only select cell lines from RNAseqcountsGlioma that are used for each group 
  dgListGlioma<- DGEList(counts=RNAseqcountsGlioma %>% 
                           select(colnames(RNAseqcountsGlioma)[3:67]==ListofAllCellLines[DesignMatrixIndex]), 
                         genes=RNAseqcountsGlioma[,1:2])
  #SD comment: This looks generally good. You may want to not hard-code the 67 on line 212 - what if one time you have
  #more cell lines than that? It would cause a problem, but there's a way to fix it.
  #Also, you may want to have your dgListGlioma to actually be a list of dgLists, one for each condition you're testing
  #That way, you'll save all the data rather than losing it each time the loop re-sets. Definitely a good start so far!
  countsPerMillion <- cpm(dgListGlioma) 
}

#Past work
# Creating a DGEList object
dgListGliomaTMZMGMTlow <- DGEList(counts=RNAseqcountsGlioma %>% select(TemozolomideMGMTlow$celllines), genes=RNAseqcountsGlioma[,1:2])
countsPerMillion <- cpm(dgListGliomaTMZMGMTlow) 
# Filtering
#cpm is number of reads mapped to a gene/total number of mapped reads from th library
countCheck <- countsPerMillion > 1 #which genes have more than 1 cpm
keep <- which(rowSums(countCheck) >= 2) #rowSums adds TRUEs for each gene, #genes to keep
dgListGliomaTMZMGMTlow <- dgListGliomaTMZMGMTlow[keep,]
# Normalization
#TMM normalization adjusts library sizes based on the assumption that most genes are not differentially expressed
dgListGliomaTMZMGMTlow <- calcNormFactors(dgListGliomaTMZMGMTlow, method="TMM") #Trimmed mean of M-values normalization
# Data Exploration
plotMDS(dgListGlioma) #what does this show again?
# Setting up the model- try TMZ/MGMT low first
designMatTemozolomideMGMTlow
# Estimating Dispersons
dgListGliomaTMZMGMTlow <- estimateGLMCommonDisp(dgListGliomaTMZMGMTlow, design=designMatTemozolomideMGMTlow)
#dispersion is same for each gene
dgListGliomaTMZMGMTlow <- estimateGLMTrendedDisp(dgListGliomaTMZMGMTlow, design=designMatTemozolomideMGMTlow)
#trend between expression and variation
dgListGliomaTMZMGMTlow <- estimateGLMTagwiseDisp(dgListGliomaTMZMGMTlow, design=designMatTemozolomideMGMTlow)
#each gene has its own dispersion
plotBCV(dgListGliomaTMZMGMTlow) 
# Differential Expression
fit <- glmFit(dgListGliomaTMZMGMTlow, designMatTemozolomideMGMTlow)
contrast_dgListGliomaTMZMGMTlow <- makeContrasts(TMZMGMTlow=resistant-sensitive,
                                                 levels=designMatTemozolomideMGMTlow)
lrt <- glmLRT(fit, contrast=contrast_dgListGliomaTMZMGMTlow)
edgeR_result <- topTags(lrt)
save(topTags(lrt,n=15000)$table, file='/Users/caitlinhalligan/Desktop/Thesis')
deGenes <- decideTestsDGE(lrt, p=0.001)
deGenes <- rownames(lrt)[as.logical(deGenes)]
plotSmear(lrt, de.tags=deGenes)
     abline(h=c(-1, 1), col=2)
#differentially expressed genes
View(edgeR_result$table)
