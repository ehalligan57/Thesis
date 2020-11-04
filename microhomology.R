#Simona Dalin
#6/23/20
#Look at relationship between microhomology length at SVs, sensitivity to alkylating agents, and MGMT status
library(tidyr)
#setwd("/Volumes/xchip_beroukhimlab/Simona/CCNU")
dirName <- getwd()
############# Load in the data ###################
#Microhomology lengths
SvABAFileNames <- list.files(sprintf("%s/ccle.svaba",dirName))
SvABAFilePaths <- file.path(dirName,"ccle.svaba",fileNames)
SvABAData <- lapply(SvABAFilePaths, readRDS)
accessedSvABAData <- sapply(SvABAData, values)
names(accessedSvABAData) <- substring(SvABAFileNames,8,nchar(SvABAFileNames)-32)
#MGMT expression
RNAseq <- fread(file.path(dirName,"ccle.RNAseq.rsem.gene","CCLE_RNAseq_rsem_genes_tpm_20180929.txt"))
######## CCLE drug data ########
CCLEGNFSensitivity <- fread(file.path(dirName,"ccle.gdsc.drugSensitivity","CCLE_GNF_data_090613.csv"))
CCLEGNFMetadata <- fread(file.path(dirName,"ccle.gdsc.drugSensitivity","CCLE_GNF_metadata_090613.csv"))
CCLENP24Metadata <- fread(file.path(dirName,"ccle.gdsc.drugSensitivity","CCLE_NP24.2009_profiling_2012.02.20.csv"))
CCLENP24Sensitivity <- fread(file.path(dirName,"ccle.gdsc.drugSensitivity","CCLE_NP24.2009_Drug_data_2015.02.24.csv"))
#No alkylators here
######## Sanger drug data ########
sangerDoseResponse <- fread(file.path(dirName,"ccle.gdsc.drugSensitivity","sanger-dose-response.csv"))
sangerViability <- fread(file.path(dirName,"ccle.gdsc.drugSensitivity","sanger-viability.csv"))
gdsc1Metadata <- fread(file.path(dirName,"ccle.gdsc.drugSensitivity","gdsc1drugmetadata.csv"))
gdsc2Metadata <- fread(file.path(dirName,"ccle.gdsc.drugSensitivity","gdsc2drugmetadata.csv"))
alkylatorRowsgdsc1 <- grep("alkyla", gdsc1Metadata$targets) #TMZ
alkylatorRowsgdsc2 <- grep("alkyla", gdsc2Metadata$targets) #"Oxaliplatin" "Oxaliplatin" "Temozolomide"
######## CTD2 drug data ##########
#("chlorambucil" "cyclophosphamide" "dacarbazine" "ifosfamide" "temozolomide" "bendamustine" "Platin" "oxaliplatin")
#Will need to look in DrugMetaData for master_cpd_id, in cellLineMetadata for master_ccl_id, and in ExperimentMetadata to match master_ccl_id to
experiment_id
CTD2DrugMetaData <- fread(file.path(dirName,"ccle.gdsc.drugSensitivity","v20.meta.per_compound.txt"))
alkylatorRowsCTD2 <- grep("alkyla", CTD2DrugMetaData$target_or_activity_of_compound)
alkylatorMaster_cpd_id <- CTD2DrugMetaData$master_cpd_id[alkylatorRowsCTD2]
CTD2CellLineMetadata <- fread(file.path(dirName,"ccle.gdsc.drugSensitivity","v20.meta.per_cell_line.txt"))
CTD2ExperimentMetadata <- fread(file.path(dirName,"ccle.gdsc.drugSensitivity","v20.meta.per_experiment.txt"))
CTD2DrugSensitivity <- fread(file.path(dirName, "ccle.gdsc.drugSensitivity", "v20.data.curves_post_qc.txt"))
############### Plot histogram of MGMT expression #################
MGMTgeneID <- "ENSG00000170430"
RNAseqMGMTrow <- grep(MGMTgeneID, RNAseq$gene_id)
MGMTexpressionValues <- RNAseq[RNAseqMGMTrow, 3:dim(RNAseq)[2]]
MGMTexpressionValues.df <- data.frame(as.numeric(MGMTexpressionValues))
MGMTexpressionValues.df$geneID <- 1

ggplot(MGMTexpressionValues.df, aes(x = geneID, y = as.numeric.MGMTexpressionValues.)) +
  geom_violin() +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.2, binwidth = 2) +
  labs(title = "MGMT expression across CCLE cell lines") +
  xlab("MGMT") +
  ylab("RSEM")
########## Extract line names of MGMT high and low cell lines ###########
highCutoff <- 19
lowCutoff <- 6
MGMTexpressionValues.df$cellLineIDs <- names(MGMTexpressionValues)
MGMTLowLines <- MGMTexpressionValues.df$cellLineNames[MGMTexpressionValues <= lowCutoff]
MGMTHighLines <- MGMTexpressionValues.df$cellLineNames[MGMTexpressionValues >= highCutoff]
#Remove the tissue annotation from cell line names
MGMTexpressionValues.df$cellLineNames <- lapply(MGMTexpressionValues.df$cellLineIDs, function(x){
  substring(x, 1, regexpr("_", x) - 1)})
########## Make Data Table with cell line, drug, MGMT expression, and drug sensitivity information ##############
#extract only drug sensitivity of compounds we're interested in
CTD2AlkylatorSensitivity <- CTD2DrugSensitivity[CTD2DrugSensitivity$master_cpd_id %in% alkylatorMaster_cpd_id,]
#Store the relevant information into a data.table
CTD2AlkylatorSensitivity.dt <- data.table(experiment_id = CTD2AlkylatorSensitivity$experiment_id,
                                          apparent_ec50_umol = CTD2AlkylatorSensitivity$apparent_ec50_umol,
                                          area_under_curve = CTD2AlkylatorSensitivity$area_under_curve,
                                          master_cpd_id = CTD2AlkylatorSensitivity$master_cpd_id)
#Add in real drug names because this cpd.id stuff is bullshit
CTD2AlkylatorSensitivity.dt$drug_name <- sapply(CTD2AlkylatorSensitivity.dt$master_cpd_id, function(x){
  CTD2DrugMetaData$cpd_name[which(x == CTD2DrugMetaData$master_cpd_id)]})
#Now translate experiment_id into cell line name...
CTD2AlkylatorSensitivity.dt$master_ccl_id <- sapply(CTD2AlkylatorSensitivity.dt$experiment_id, function(x){
  CTD2ExperimentMetadata$master_ccl_id[which(x == CTD2ExperimentMetadata$experiment_id)]})
CTD2AlkylatorSensitivity.dt$ccl_name <- sapply(CTD2AlkylatorSensitivity.dt$master_ccl_id, function(x){
  CTD2CellLineMetadata$ccl_name[which(x == CTD2CellLineMetadata$master_ccl_id)]})
#Add MGMT expression to the data table
CTD2AlkylatorSensitivity.dt$MGMT_expression <- as.numeric(sapply(CTD2AlkylatorSensitivity.dt$ccl_name, function(x){
  MGMTexpressionValues.df$as.numeric.MGMTexpressionValues.[which(x == MGMTexpressionValues.df$cellLineNames)]}))
#Add MGMT high/low label
CTD2AlkylatorSensitivity.dt$MGMT_hi.lo[CTD2AlkylatorSensitivity.dt$MGMT_expression <= lowCutoff] <- "low"
CTD2AlkylatorSensitivity.dt$MGMT_hi.lo[CTD2AlkylatorSensitivity.dt$MGMT_expression >= highCutoff] <- "high"
CTD2AlkylatorSensitivity.dt$MGMT_hi.lo[is.na(CTD2AlkylatorSensitivity.dt$MGMT_hi.lo)] <- "moderate"
########### Plot MGMT expression vs. alkylator sensitivity ##############
#AUC
linearModels <- unique(CTD2AlkylatorSensitivity.dt$drug_name)
linearModels <- data.table(drugNames = linearModels,
                           models = NA)
for(drug in linearModels$drugNames){
  sensitivityThisDrug <- CTD2AlkylatorSensitivity.dt$area_under_curve[CTD2AlkylatorSensitivity.dt$drug_name == drug]
  MGMTThisDrug <- CTD2AlkylatorSensitivity.dt$MGMT_expression[CTD2AlkylatorSensitivity.dt$drug_name == drug]
  linearModels$models[linearModels$drugNames == drug] <- list(lm(sensitivityThisDrug ~ MGMTThisDrug))
}
MGMTexpressionAlkylatorSensitivity <- ggplot(data = CTD2AlkylatorSensitivity.dt,
                                             aes(x = MGMT_expression, y = area_under_curve, color = drug_name, shape = MGMT_hi.lo)) +
  geom_point(size = 1) +
  scale_color_discrete() +
  geom_smooth(method = lm) +
  labs(title = "MGMT expression vs. alkylator sensitivity, AUC")
ggsave("alkylatorSensitivityVsMGMTexpressionAUC.pdf",MGMTexpressionAlkylatorSensitivity, device = "pdf",height = 10, width = 10)
#EC50
linearModels <- unique(CTD2AlkylatorSensitivity.dt$drug_name)
linearModels <- data.table(drugNames = linearModels,
                           models = NA)
for(drug in linearModels$drugNames){
  sensitivityThisDrug <- CTD2AlkylatorSensitivity.dt$apparent_ec50_umol[CTD2AlkylatorSensitivity.dt$drug_name == drug]
  MGMTThisDrug <- CTD2AlkylatorSensitivity.dt$MGMT_expression[CTD2AlkylatorSensitivity.dt$drug_name == drug]
  linearModels$models[linearModels$drugNames == drug] <- list(lm(sensitivityThisDrug ~ MGMTThisDrug))
}
CTD2AlkylatorSensitivity.dt$apparent_ec50_umol[CTD2AlkylatorSensitivity.dt$apparent_ec50_umol > 200] <- 200
MGMTexpressionAlkylatorSensitivity <- ggplot(data = CTD2AlkylatorSensitivity.dt,
                                             aes(x = MGMT_expression, y = apparent_ec50_umol, color = drug_name, shape = MGMT_hi.lo)) +
  geom_point(size = 1) +
  scale_color_discrete() +
  geom_smooth(method = lm) +
  labs(title = "MGMT expression vs. alkylator sensitivity, EC50")
ggsave("alkylatorSensitivityVsMGMTexpressionEC50.pdf",MGMTexpressionAlkylatorSensitivity, device = "pdf",height = 10, width = 10)
######### Make histogram of AUC and ec50 for each drug and call each cell line as sensitive or resistant #############
#Assign each cell line as sensitive or resistant to a particular drug
CTD2AlkylatorSensitivity.dt$alkylatorSensitivityClass <- NA
for(drug in unique(CTD2AlkylatorSensitivity.dt$drug_name)){
  
  thisDrugIndices <- CTD2AlkylatorSensitivity.dt$drug_name == drug
  
  quartiles <- summary(CTD2AlkylatorSensitivity.dt$area_under_curve[thisDrugIndices])
  CTD2AlkylatorSensitivity.dt$alkylatorSensitivityClass[CTD2AlkylatorSensitivity.dt$area_under_curve <= quartiles["1st Qu."] & thisDrugIndices] <-
    "sensitive"
  CTD2AlkylatorSensitivity.dt$alkylatorSensitivityClass[CTD2AlkylatorSensitivity.dt$area_under_curve >= quartiles["3rd Qu."] & thisDrugIndices] <-
    "resistant"
  CTD2AlkylatorSensitivity.dt$alkylatorSensitivityClass[is.na(CTD2AlkylatorSensitivity.dt$alkylatorSensitivityClass) & thisDrugIndices] <- "moderate"
}
#AUCs
AUCViolinPlots <- ggplot(CTD2AlkylatorSensitivity.dt, aes(x = drug_name, y = area_under_curve)) +
  geom_violin() +
  geom_dotplot(aes(x = drug_name, y = area_under_curve, color = alkylatorSensitivityClass),
               binaxis = "y", stackdir = "center", dotsize = 0.02, binwidth = 0.2) +
  labs(x = "Drug",
       y = "AUC",
       title = "Alkylating agent sensitivity") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 10)) +
  theme(axis.title.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) +
  theme(axis.text.y = element_text(size = 10))
ggsave("AUCviolinPlots.pdf",AUCViolinPlots, device = "pdf",height = 5, width = 5)
#EC50s
EC50ViolinPlots <- ggplot(CTD2AlkylatorSensitivity.dt, aes(x = drug_name, y = apparent_ec50_umol)) +
  geom_violin() +
  geom_dotplot(aes(x = drug_name, y = area_under_curve, color = alkylatorSensitivityClass),
               binaxis = "y", stackdir = "center", dotsize = 0.2, binwidth = 0.2) +
  labs(x = "Drug",
       y = "Apparent EC50 (umol)",
       title = "Alkylating agent sensitivity") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 10)) +
  theme(axis.title.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) +
  theme(axis.text.y = element_text(size = 10))
ggsave("EC50violinPlots.pdf",EC50ViolinPlots, device = "pdf",height = 5, width = 5)
######### Extract microhomology distributions for each cell line ##############
#assign homology lengths, since they're missing
for(cellLine in names(accessedSvABAData)){
  accessedSvABAData[[cellLine]]$HOMLEN <- nchar(accessedSvABAData[[cellLine]]$HOMSEQ)
  
  alkylatorDataTableRowsThisCellLine <- CTD2AlkylatorSensitivity.dt$ccl_name == cellLine
  CTD2AlkylatorSensitivity.dt$microhomologyLengths[alkylatorDataTableRowsThisCellLine] <- list(accessedSvABAData[[cellLine]]$HOMLEN)
}
########### Extract microhomology distributions for res/sens lines in MGMT-low lines ###########
microhomologyHistograms <- data.table(drugs = NA,
                                      sensitivityClass = NA,
                                      bins = NA,
                                      counts = NA,
                                      MGMTstatus = NA,
                                      density = NA)
MGMTlowData <- CTD2AlkylatorSensitivity.dt[MGMT_hi.lo == "low"]
tempMicrohomologyHistogramResult <- plotMicrohomologyHistograms(MGMTlowData, MGMTstatus = "low")
microhomologyHistograms <- rbind(microhomologyHistograms, tempMicrohomologyHistogramResult)
MGMTmoderateData <- CTD2AlkylatorSensitivity.dt[MGMT_hi.lo == "moderate"]
tempMicrohomologyHistogramResult <- plotMicrohomologyHistograms(MGMTlowData, MGMTstatus = "moderate")
microhomologyHistograms <- rbind(microhomologyHistograms, tempMicrohomologyHistogramResult)
MGMThighData <- CTD2AlkylatorSensitivity.dt[MGMT_hi.lo == "high"]
tempMicrohomologyHistogramResult <- plotMicrohomologyHistograms(MGMTlowData, MGMTstatus = "high")
microhomologyHistograms <- rbind(microhomologyHistograms, tempMicrohomologyHistogramResult)
plotMicrohomologyHistograms <- function(MGMTdata, MGMTstatus){
  
  allMicrohomologyHistograms <- data.table(drugs = NA,
                                           sensitivityClass = NA,
                                           bins = NA,
                                           counts = NA,
                                           MGMTstatus = NA,
                                           density = NA)
  
  microhomologyDistributions <- data.table(drugs = NA,
                                           sensitivityClass = NA,
                                           microhomologyLength = NA,
                                           MGMTStatus = NA)
  
  for(drug in unique(MGMTdata$drug_name)){
    #sensitiveMicrohomologies <- as.numeric(unlist(MGMTdata[drug_name == drug][alkylatorSensitivityClass == "sensitive"]))
    #resistantMicrohomologies <- as.numeric(unlist(MGMTdata[drug_name == drug][alkylatorSensitivityClass == "resistant"]))
    #moderateMicrohomologies <- as.numeric(unlist(MGMTdata[drug_name == drug][alkylatorSensitivityClass == "moderate"]))
    
    
    for(sensitivityClass in c("sensitive", "resistant", "moderate")){
      
      if(exists("microhomologyLengthsThisClass")){rm(microhomologyLengthsThisClass)}
      microhomologyLengthsThisClass <- unlist(MGMTdata[drug_name == drug][alkylatorSensitivityClass == sensitivityClass]$microhomologyLengths)
      microhomologyLengthsThisClass[is.na(microhomologyLengthsThisClass)] <- 0
      
      microhomologyDistributionsTemp <- data.table(drugs = drug,
                                                   sensitivityClass = sensitivityClass,
                                                   microhomologyLength = microhomologyLengthsThisClass,
                                                   MGMTStatus = MGMTstatus)
      
      microhomologyDistributions <- rbind(microhomologyDistributions, microhomologyDistributionsTemp)
      
    }
    
    
    #Remove microhomology lengths greater than 51bp
    microhomologyDistributionsLessThan51bp <- microhomologyDistributions[microhomologyDistributions$microhomologyLength < 51,]
    
    #Make histograms for each type of drug sensitivity
    sensitiveHistogram <- hist(microhomologyDistributionsLessThan51bp$microhomologyLength[microhomologyDistributionsLessThan51bp$sensitivityClass ==
                                                                                            "sensitive"],
                               breaks = c(-1:51))
    sensitive.dt <- data.table(drugs = drug,
                               sensitivityClass = "sensitive",
                               bins = sensitiveHistogram$breaks[1:51]+1,
                               counts = sensitiveHistogram$counts[1:51],
                               density = sensitiveHistogram$density[1:51],
                               MGMTstatus = MGMTstatus)
    
    resistantHistogram <- hist(microhomologyDistributionsLessThan51bp$microhomologyLength[microhomologyDistributionsLessThan51bp$sensitivityClass ==
                                                                                            "resistant"],
                               breaks = c(-1:51))
    resistant.dt <- data.table(drugs = drug,
                               sensitivityClass = "resistant",
                               bins = resistantHistogram$breaks[1:51]+1,
                               counts = resistantHistogram$counts[1:51],
                               density = resistantHistogram$density[1:51],
                               MGMTstatus = MGMTstatus)
    
    moderateHistogram <- hist(microhomologyDistributionsLessThan51bp$microhomologyLength[microhomologyDistributionsLessThan51bp$sensitivityClass ==
                                                                                           "moderate"],
                              breaks = c(-1:51))
    moderate.dt <- data.table(drugs = drug,
                              sensitivityClass = "moderate",
                              bins = moderateHistogram$breaks[1:51]+1,
                              counts = moderateHistogram$counts[1:51],
                              density = moderateHistogram$density[1:51],
                              MGMTstatus = MGMTstatus)
    
    microhomologyHistograms <- rbind(sensitive.dt, resistant.dt, moderate.dt)
    allMicrohomologyHistograms <- rbind(allMicrohomologyHistograms, microhomologyHistograms)
    
    
    #Plot histograms of microhomology distributions for each drug by count
    counts <- ggplot(microhomologyHistograms[microhomologyHistograms$drugs == drug,], aes(x = bins, y = counts, group = sensitivityClass)) +
      geom_line(aes(color = sensitivityClass)) +
      geom_point(aes(color = sensitivityClass)) +
      labs(x = "bases of microhomology",
           y = "Count of microhomology bases",
           title = sprintf("Microhomology distribution in MGMT %s lines\n%s", MGMTstatus, drug)) +
      theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
      theme(axis.title.x = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(axis.text.x = element_text(size = 10)) +
      theme(axis.text.y = element_text(size = 10))
    
    
    ggsave(sprintf("%s Microhomology Distributions in MGMT %s Lines By Count.pdf", drug, MGMTstatus),counts, device = "pdf",height = 5, width = 10)
    
    
    #Plot histograms of microhomology distributions for each drug by density
    density <- ggplot(microhomologyHistograms[microhomologyHistograms$drugs == drug,], aes(x = bins, y = density, group = sensitivityClass)) +
      geom_line(aes(color = sensitivityClass)) +
      geom_point(aes(color = sensitivityClass)) +
      labs(x = "bases of microhomology",
           y = "Density of microhomology bases",
           title = sprintf("Microhomology distribution in MGMT %s lines\n%s", MGMTstatus, drug)) +
      theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
      theme(axis.title.x = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(axis.text.x = element_text(size = 10)) +
      theme(axis.text.y = element_text(size = 10))
    
    
    ggsave(sprintf("%s Microhomology Distributions in MGMT %s Lines By Density.pdf", drug, MGMTstatus),density, device = "pdf",height = 5, width = 10)
    
  }
  
  return(allMicrohomologyHistograms)
}
#Remove NA rows from microhomology histograms...
microhomologyHistograms <- microhomologyHistograms[!apply(microhomologyHistograms, MARGIN = 1, function(x){all(is.na(x))})]
save(CTD2AlkylatorSensitivity.dt, file = "CTD2AlkylatorSensitivity.Rdata")
save(microhomologyHistograms, file = "microhomologyHistograms.Rdata")

########### Subset microhomology distributions into DSB repair types and perform statistical test ###########
temp <- microhomologyHistograms[,.N, by=.(drugs,MGMTstatus)][,1:2]
microhomologyBins <- data.table(temp, Counts = NA)
for(drugMGMTstatus in 1:dim(microhomologyBins)[1]){
  drug <- microhomologyBins$drugs[drugMGMTstatus]
  MGMTstatus <- microhomologyBins$MGMTstatus[drugMGMTstatus]
  
  #Extract histogram counts for this drug/mgmtstatus
  rowsToExtract <- microhomologyHistograms$drugs == drug & microhomologyHistograms$MGMTstatus == MGMTstatus
  tempHistogram <- microhomologyHistograms[rowsToExtract,]
  setkey(tempHistogram, bins, sensitivityClass)
  
  binnedDSBrepair <- data.table(matrix(nrow = 3, ncol = 2))
  colnames(binnedDSBrepair) <- c("sensitive","resistant")
  binnedDSBrepair <- apply(binnedDSBrepair, MARGIN = 2, function(x){as.integer(x)})
  rownames(binnedDSBrepair) <- c("NHEJ","MMEJ","SSA")
  
  
  for(sensitivityClassName in c("sensitive","resistant")){
    sensitivityClassCol <- colnames(binnedDSBrepair) == sensitivityClassName
    
    #NHEJ
    bases <- 0
    binnedDSBrepair[1, sensitivityClassCol] <- tempHistogram[.(bases,sensitivityClassName), sum(counts)]
    
    #MMEJ
    bases <- 4:9
    binnedDSBrepair[2, sensitivityClassCol] <- tempHistogram[.(bases,sensitivityClassName), sum(counts)]
    
    #SSA
    bases <- 13:50
    binnedDSBrepair[3, sensitivityClassCol] <- tempHistogram[.(bases,sensitivityClassName), sum(counts)]
  }
  
  microhomologyBins$Counts[drugMGMTstatus] <- list(binnedDSBrepair)
  
  
  chiSquaredResults <- chisq.test(binnedDSBrepair)
  microhomologyBins$chisqTest[drugMGMTstatus] <- list(chiSquaredResults)
  
  #Plot residuals if the test was significant
  if(chiSquaredResults$p.value < (0.05/nrow(microhomologyBins))){
    residuals <- as.data.frame(chiSquaredResults$residuals)
    residuals$`DSB Repair` <- rownames(residuals)
    tidyResiduals <- residuals %>%
      pivot_longer(-`DSB Repair`, names_to = "Alkylator Sensitivity", values_to = "Residuals")
    
    #Sort tidyResiduals so it shows up in the right order on the plot
    tidyResiduals$`DSB Repair` <- factor(tidyResiduals$`DSB Repair`, levels = c("NHEJ","MMEJ","SSA"))
    
    #Plot
    residualsPlot <- ggplot(tidyResiduals, aes(x = `DSB Repair`, y = Residuals, fill = `Alkylator Sensitivity`)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      labs(title = sprintf("Residuals of chi-squared comparing microhomology distributions in \n%s sensitive and resistant MGMT %s CCLE lines", drug,
                           MGMTstatus),
           subtitle = sprintf("p-value = %.2f", chiSquaredResults$p.value)) +
      theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
      theme(axis.title.x = element_text(size = 15)) +
      theme(axis.title.y = element_text(size = 15)) +
      theme(axis.text.x = element_text(size = 15)) +
      theme(axis.text.y = element_text(size = 15))
    
    
    ggsave(sprintf("residualsTo%sinMGMT%sLines.pdf",drug, MGMTstatus), residualsPlot, device = "pdf",height = 10, width = 10)