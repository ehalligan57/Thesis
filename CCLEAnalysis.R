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
for(ColNumber in (3:ncol(EdgeRdataframe))){#loops over drugs
  MGMThighSensitivity <- c()
  MGMTlowSensitivity <- c()
  for(RowNumber in (1:nrow(EdgeRdataframe))){ #loops of cell lines
    if (EdgeRdataframe$MGMTexpressiondata[RowNumber]=='high'){
      if (is.na(EdgeRdataframe[RowNumber,ColNumber])){ #skip NAs
        next #This is the "official" command to skip to the next iteration of a for or while loop. It does work without it too.
      }else{
        new_sensitivity <- EdgeRdataframe[RowNumber,ColNumber]
        MGMThighSensitivity <- c(MGMThighSensitivity, new_sensitivity)
      }
    }else{ #MGMT low
      if (is.na(EdgeRdataframe[RowNumber,ColNumber])){ #skip NAs
      }else{
        new_sensitivity <- EdgeRdataframe[RowNumber,ColNumber]
        MGMTlowSensitivity <- c(MGMTlowSensitivity, new_sensitivity)
      }
    }
  }
  new_designMatMGMThigh <- model.matrix(~0+MGMThighSensitivity) #this works, same as old one
  new_designMatMGMTlow <- model.matrix(~0+MGMTlowSensitivity) #this works, same as old one
  
  #saying "ListofDesignMatricesMGMThigh[[ColNumber]]" was storing each
  #design matrix into the nth + 2 slot, since the first colNumber
  #you use was col #3! I changed it so they start at 1 now.
  nextEmptySlot <- length(ListofDesignMatricesMGMThigh) + 1
  
  ListofDesignMatricesMGMThigh[[nextEmptySlot]] <- new_designMatMGMThigh
  ListofDesignMatricesMGMTlow[[nextEmptySlot]] <- new_designMatMGMTlow
} 

#Name the slots of the lists of design matrices by the drug they refer to
names(ListofDesignMatricesMGMThigh) <- as.vector(sapply(names(EdgeRdataframe)[3:ncol(EdgeRdataframe)], function(x){substr(x, 1, regexpr(" ", x)-1)}))
names(ListofDesignMatricesMGMTlow) <- sapply(names(EdgeRdataframe)[3:ncol(EdgeRdataframe)], function(x){substr(x, 1, regexpr(" ", x)-1)})



ListofDesignMatricesMGMThigh[[3]]==designMatChlorambucilMGMThigh 
#it's the same for Chlorambucil, Dacarbazine, Temozolomide, Bendamustine, Platin, Cyclophosphamide, Oxaliplatin
#it is wrong for Ifosfamide (5)

#this is true for me:
all(designMatIfosfamideMGMThigh == ListofDesignMatricesMGMThigh$ifosfamide)

ListofDesignMatricesMGMTlow[[3]]==designMatChlorambucilMGMTlow
#it's right for Dacarbazine, Temozolomide, Bendamustine, Platin, Cyclophosphamide, Oxaliplatin
#it's wrong for Chlorambucil (3), Ifosfamide (5)

#this is true for me:
all(designMatChlorambucilMGMTlow == ListofDesignMatricesMGMTlow$chlorambucil)

#This gives me an error:
all(designMatIfosfamideMGMTlow == ListofDesignMatricesMGMTlow$ifosfamide)

#They are different sizes:
dim(designMatIfosfamideMGMTlow)
dim(ListofDesignMatricesMGMTlow$ifosfamide)

#There are 14 MGMT low cell lines with ifosphamide sensitivity data
sum(!is.na(EdgeRdataframe$`ifosfamide SensitivityAUC`)[EdgeRdataframe$MGMTexpressiondata == 'low'])

#So the new version with the for loop is correct. Why was the old version incorrect though?
#There's a typo on line 290 where you set that design matrix referring to MGMT high, not low.


# Design Matrix for MGMT high, chlorambucil # #28 total, 14 R, 14S #error message but it looks right?
ChlorambucilMGMThigh <- EdgeRdataframe %>%
  dplyr::filter(MGMTexpressiondata=='high') %>%
  select(celllines,`chlorambucil SensitivityAUC`) %>%
  na.omit(`chlorambucil SensitivityAUC`)
designMatChlorambucilMGMThigh <- model.matrix(~0+ChlorambucilMGMThigh$`chlorambucil SensitivityAUC`)

# Design Matrix for MGMT low, chlorambucil #15 total, 8 R, 7 S
#ChlorambucilMGMTlow <- EdgeRdataframe %>% filter(EdgeRdataframe$MGMTexpression=='low')
#ChlorambucilMGMTlow <- data.frame(celllines=ChlorambucilMGMTlow$celllines,ChlorambucilMGMTlow = ChlorambucilMGMTlow$`chlorambucil SensitivityAUC`)
#ChlorambucilMGMTlow <- na.omit(ChlorambucilMGMTlow)
#designMatChlorambucilMGMTlow <- model.matrix(~ChlorambucilMGMTlow$ChlorambucilMGMTlow)
ChlorambucilMGMTlow <- EdgeRdataframe %>%
  dplyr::filter(MGMTexpressiondata=='low') %>%
  select(celllines,`chlorambucil SensitivityAUC`) %>%
  na.omit(`chlorambucil SensitivityAUC`)
designMatChlorambucilMGMTlow <- model.matrix(~0+ChlorambucilMGMTlow$`chlorambucil SensitivityAUC`)

# Design Matrix for MGMT high, dacarbazine # #28 total, 12 R, 16 S
#DacarbazineMGMThigh <- EdgeRdataframe %>% filter(EdgeRdataframe$MGMTexpression=='high')
#DacarbazineMGMThigh <- data.frame(celllines=DacarbazineMGMThigh$celllines,DacarbazineMGMThigh = DacarbazineMGMThigh$`dacarbazine SensitivityAUC`)
#DacarbazineMGMThigh <- na.omit(DacarbazineMGMThigh)
#designMatDacarbazineMGMThigh <- model.matrix(~DacarbazineMGMThigh$DacarbazineMGMThigh)
DacarbazineMGMThigh <- EdgeRdataframe %>%
  dplyr::filter(MGMTexpressiondata=='high') %>%
  select(celllines,`dacarbazine SensitivityAUC`) %>%
  na.omit(`dacarbazine SensitivityAUC`)
designMatDacarbazineMGMThigh <- model.matrix(~0+DacarbazineMGMThigh$`dacarbazine SensitivityAUC`)

# Design Matrix for MGMT low, dacarbazine # #15 total, 10 R, 5 S
#DacarbazineMGMTlow <- EdgeRdataframe %>% filter(EdgeRdataframe$MGMTexpression=='low')
#DacarbazineMGMTlow <- data.frame(celllines=DacarbazineMGMTlow$celllines,DacarbazineMGMTlow = DacarbazineMGMTlow$`dacarbazine SensitivityAUC`)
#DacarbazineMGMTlow <- na.omit(DacarbazineMGMTlow)
#designMatDacarbazineMGMTlow <- model.matrix(~DacarbazineMGMTlow$DacarbazineMGMTlow)
DacarbazineMGMTlow <- EdgeRdataframe %>%
  dplyr::filter(MGMTexpressiondata=='low') %>%
  select(celllines,`dacarbazine SensitivityAUC`) %>%
  na.omit(`dacarbazine SensitivityAUC`)
designMatDacarbazineMGMTlow <- model.matrix(~0+DacarbazineMGMTlow$`dacarbazine SensitivityAUC`)

# Design Matrix for MGMT high, ifosfamide # #13 total, 6 R, 7 S
#IfosfamideMGMThigh <- EdgeRdataframe %>% filter(EdgeRdataframe$MGMTexpression=='high')
#IfosfamideMGMThigh <- data.frame(celllines=IfosfamideMGMThigh$celllines,IfosfamideMGMThigh = IfosfamideMGMThigh$`ifosfamide SensitivityAUC`)
#IfosfamideMGMThigh <- na.omit(IfosfamideMGMThigh)
#designMatIfosfamideMGMThigh <- model.matrix(~IfosfamideMGMThigh$IfosfamideMGMThigh)

IfosfamideMGMThigh <- EdgeRdataframe %>%
  dplyr::filter(MGMTexpressiondata=='high') %>%
  select(celllines,`ifosfamide SensitivityAUC`) %>%
  na.omit(`ifosfamide SensitivityAUC`)
designMatIfosfamideMGMThigh <- model.matrix(~0+IfosfamideMGMThigh$`ifosfamide SensitivityAUC`)

# Design Matrix for MGMT low, ifosfamide # #9 total, 5 R, 4 S
#IfosfamideMGMTlow <- EdgeRdataframe %>% filter(EdgeRdataframe$MGMTexpression=='low')
#IfosfamideMGMTlow <- data.frame(celllines=IfosfamideMGMTlow$celllines,IfosfamideMGMTlow = IfosfamideMGMTlow$`ifosfamide SensitivityAUC`)
#IfosfamideMGMTlow <- na.omit(IfosfamideMGMTlow)
#designMatIfosfamideMGMTlow <- model.matrix(~IfosfamideMGMTlow$IfosfamideMGMTlow)

IfosfamideMGMTlow <- EdgeRdataframe %>%
  dplyr::filter(MGMTexpressiondata=='low') %>%
  select(celllines,`ifosfamide SensitivityAUC`) %>%
  na.omit(`ifosfamide SensitivityAUC`)
designMatIfosfamideMGMTlow <- model.matrix(~0+IfosfamideMGMThigh$`ifosfamide SensitivityAUC`) #typo here: referring to Ifosfamide high, not low.

# Design Matrix for MGMT high, temozolomide # #28 total, 17 R, 11 S
#TemozolomideMGMThigh <- EdgeRdataframe %>% filter(EdgeRdataframe$MGMTexpression=='high')
#TemozolomideMGMThigh <- data.frame(celllines=TemozolomideMGMThigh$celllines,TemozolomideMGMThigh = TemozolomideMGMThigh$`temozolomide SensitivityAUC`)
#TemozolomideMGMThigh <- na.omit(TemozolomideMGMThigh)
#designMatTemozolomideMGMThigh <- model.matrix(~TemozolomideMGMThigh$TemozolomideMGMThigh)

TemozolomideMGMThigh <- EdgeRdataframe %>%
  dplyr::filter(MGMTexpressiondata=='high') %>%
  select(celllines,`temozolomide SensitivityAUC`) %>%
  na.omit(`temozolomide SensitivityAUC`)
designMatTemozolomideMGMThigh <- model.matrix(~0+TemozolomideMGMThigh$`temozolomide SensitivityAUC`)

# Design Matrix for MGMT low, temozolomide # #17 total, 9 R, 8 S
#TemozolomideMGMTlow <- EdgeRdataframe %>% filter(EdgeRdataframe$MGMTexpression=='low')
#TemozolomideMGMTlow <- data.frame(celllines=TemozolomideMGMTlow$celllines,TemozolomideMGMTlow = TemozolomideMGMTlow$`temozolomide SensitivityAUC`)
#TemozolomideMGMTlow <- na.omit(TemozolomideMGMTlow)
#designMatTemozolomideMGMTlow <- model.matrix(~TemozolomideMGMTlow$TemozolomideMGMTlow)
TemozolomideMGMTlow <- EdgeRdataframe %>%
  dplyr::filter(MGMTexpressiondata=='low') %>%
  select(celllines,`temozolomide SensitivityAUC`) %>%
  na.omit(`temozolomide SensitivityAUC`)
designMatTemozolomideMGMTlow <- model.matrix(~0+TemozolomideMGMTlow$`temozolomide SensitivityAUC`)
colnames(designMatTemozolomideMGMTlow) <- c("resistant","sensitive")
# Design Matrix for MGMT high, bendamustine # #27 total, 14 R, 13 S
#BendamustineMGMThigh <- EdgeRdataframe %>% filter(EdgeRdataframe$MGMTexpression=='high')
#BendamustineMGMThigh <- data.frame(celllines=BendamustineMGMThigh$celllines,BendamustineMGMThigh = BendamustineMGMThigh$`bendamustine SensitivityAUC`)
#BendamustineMGMThigh <- na.omit(BendamustineMGMThigh)
#designMatBendamustineMGMThigh <- model.matrix(~BendamustineMGMThigh$BendamustineMGMThigh)
BendamustineMGMThigh <- EdgeRdataframe %>%
  dplyr::filter(MGMTexpressiondata=='high') %>%
  select(celllines,`bendamustine SensitivityAUC`) %>%
  na.omit(`bendamustine SensitivityAUC`)
designMatBendamustineMGMThigh <- model.matrix(~0+BendamustineMGMThigh$`bendamustine SensitivityAUC`)

# Design Matrix for MGMT low, bendamustine # #17 total, 8 R, 9 S
#BendamustineMGMTlow <- EdgeRdataframe %>% filter(EdgeRdataframe$MGMTexpression=='low')
#BendamustineMGMTlow <- data.frame(celllines=BendamustineMGMTlow$celllines,BendamustineMGMTlow = BendamustineMGMTlow$`bendamustine SensitivityAUC`)
#BendamustineMGMTlow <- na.omit(BendamustineMGMTlow)
#designMatBendamustineMGMTlow <- model.matrix(~BendamustineMGMTlow$BendamustineMGMTlow)
BendamustineMGMTlow <- EdgeRdataframe %>%
  dplyr::filter(MGMTexpressiondata=='low') %>%
  select(celllines,`bendamustine SensitivityAUC`) %>%
  na.omit(`bendamustine SensitivityAUC`)
designMatBendamustineMGMTlow <- model.matrix(~0+BendamustineMGMTlow$`bendamustine SensitivityAUC`)

# Design Matrix for MGMT high, Platin # #29 total, 15 R, 14 S
#PlatinMGMThigh <- EdgeRdataframe %>% filter(EdgeRdataframe$MGMTexpression=='high')
#PlatinMGMThigh <- data.frame(celllines=PlatinMGMThigh$celllines,PlatinMGMThigh = PlatinMGMThigh$`Platin SensitivityAUC`)
#PlatinMGMThigh <- na.omit(PlatinMGMThigh)
#designMatPlatinMGMThigh <- model.matrix(~PlatinMGMThigh$PlatinMGMThigh)
PlatinMGMThigh <- EdgeRdataframe %>%
  dplyr::filter(MGMTexpressiondata=='high') %>%
  select(celllines,`Platin SensitivityAUC`) %>%
  na.omit(`Platin SensitivityAUC`)
designMatPlatinMGMThigh <- model.matrix(~0+PlatinMGMThigh$`Platin SensitivityAUC`)

# Design Matrix for MGMT low, Platin # #17 total, 8 R, 9 S
#PlatinMGMTlow <- EdgeRdataframe %>% filter(EdgeRdataframe$MGMTexpression=='low')
#PlatinMGMTlow <- data.frame(celllines=PlatinMGMTlow$celllines,PlatinMGMTlow = PlatinMGMTlow$`Platin SensitivityAUC`)
#PlatinMGMTlow <- na.omit(PlatinMGMTlow)
#designMatPlatinMGMTlow <- model.matrix(~PlatinMGMTlow$PlatinMGMTlow)
PlatinMGMTlow <- EdgeRdataframe %>%
  dplyr::filter(MGMTexpressiondata=='low') %>%
  select(celllines,`Platin SensitivityAUC`) %>%
  na.omit(`Platin SensitivityAUC`)
designMatPlatinMGMTlow <- model.matrix(~0+PlatinMGMTlow$`Platin SensitivityAUC`)

# Design Matrix for MGMT high, cyclophosphamide # #
#CyclophosphamideMGMThigh <- EdgeRdataframe %>% filter(EdgeRdataframe$MGMTexpression=='high')
#CyclophosphamideMGMThigh <- data.frame(celllines=CyclophosphamideMGMThigh$celllines,CyclophosphamideMGMThigh = CyclophosphamideMGMThigh$`cyclophosphamide SensitivityAUC`)
#CyclophosphamideMGMThigh <- na.omit(CyclophosphamideMGMThigh)
#designMatCyclophosphamideMGMThigh <- model.matrix(~CyclophosphamideMGMThigh$CyclophosphamideMGMThigh)
CyclophosphamideMGMThigh <- EdgeRdataframe %>%
  dplyr::filter(MGMTexpressiondata=='high') %>% 
  select(celllines,`cyclophosphamide SensitivityAUC`) %>%
  na.omit(`cyclophosphamide SensitivityAUC`)
designMatCyclophosphamideMGMThigh <- model.matrix(~0+CyclophosphamideMGMThigh$`cyclophosphamide SensitivityAUC`)

# Design Matrix for MGMT low, cyclophosphamide # #9 total, 4 R, 5 S
#CyclophosphamideMGMTlow <- EdgeRdataframe %>% filter(EdgeRdataframe$MGMTexpression=='low')
#CyclophosphamideMGMTlow <- data.frame(celllines=CyclophosphamideMGMTlow$celllines,CyclophosphamideMGMTlow = CyclophosphamideMGMTlow$`cyclophosphamide SensitivityAUC`)
#CyclophosphamideMGMTlow <- na.omit(CyclophosphamideMGMTlow)
#designMatCyclophosphamideMGMTlow <- model.matrix(~CyclophosphamideMGMTlow$CyclophosphamideMGMTlow)
CyclophosphamideMGMTlow <- EdgeRdataframe %>%
  dplyr::filter(MGMTexpressiondata=='low') %>% 
  select(celllines,`cyclophosphamide SensitivityAUC`) %>%
  na.omit(`cyclophosphamide SensitivityAUC`)
designMatCyclophosphamideMGMTlow <- model.matrix(~0+CyclophosphamideMGMTlow$`cyclophosphamide SensitivityAUC`)

# Design Matrix for MGMT high, oxaliplatin # #28 total, 14 R, 14 S
#OxaliplatinMGMThigh <- EdgeRdataframe %>% filter(EdgeRdataframe$MGMTexpression=='high')
#OxaliplatinMGMThigh <- data.frame(celllines=OxaliplatinMGMThigh$celllines,OxaliplatinMGMThigh = OxaliplatinMGMThigh$`oxaliplatin SensitivityAUC`)
#OxaliplatinMGMThigh <- na.omit(OxaliplatinMGMThigh)
#designMatOxaliplatinMGMThigh <- model.matrix(~OxaliplatinMGMThigh$OxaliplatinMGMThigh)
OxaliplatinMGMThigh <- EdgeRdataframe %>%
  dplyr::filter(MGMTexpressiondata=='high') %>% 
  select(celllines,`oxaliplatin SensitivityAUC`) %>%
  na.omit(`oxaliplatin SensitivityAUC`)
designMatOxaliplatinMGMThigh <- model.matrix(~0+OxaliplatinMGMThigh$`oxaliplatin SensitivityAUC`)

# Design Matrix for MGMT low, oxaliplatin # #17 total, 8 R, 8 S
#OxaliplatinMGMTlow <- EdgeRdataframe %>% filter(EdgeRdataframe$MGMTexpression=='low')
#OxaliplatinMGMTlow <- data.frame(celllines=OxaliplatinMGMTlow$celllines,OxaliplatinMGMTlow = OxaliplatinMGMTlow$`oxaliplatin SensitivityAUC`)
#OxaliplatinMGMTlow <- na.omit(OxaliplatinMGMTlow)
#designMatOxaliplatinMGMTlow <- model.matrix(~OxaliplatinMGMTlow$OxaliplatinMGMTlow)
OxaliplatinMGMTlow <- EdgeRdataframe %>%
  dplyr::filter(MGMTexpressiondata=='low') %>% 
  select(celllines,`oxaliplatin SensitivityAUC`) %>%
  na.omit(`oxaliplatin SensitivityAUC`)
designMatOxaliplatinMGMTlow <- model.matrix(~0+OxaliplatinMGMTlow$`oxaliplatin SensitivityAUC`)

#########     EdgeR Workflow    ######### 
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
# why do the rows and col have to be the same?
# estimates dispersions (how far apart the data is), 3 different methods
# use common dispersions, an estimate based on mean-variance, genewise dispersion
# takes into account dispersion when calculating significant differential expression
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
#contrast vs co_eff?
edgeR_result <- topTags(lrt)
save(topTags(lrt,n=15000)$table, file='/Users/caitlinhalligan/Desktop/Thesis')
deGenes <- decideTestsDGE(lrt, p=0.001)
deGenes <- rownames(lrt)[as.logical(deGenes)]
plotSmear(lrt, de.tags=deGenes)
     abline(h=c(-1, 1), col=2)
#differentially expressed genes
View(edgeR_result$table)
