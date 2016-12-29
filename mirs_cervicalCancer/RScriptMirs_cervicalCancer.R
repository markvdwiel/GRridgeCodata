# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                       Load packages                                #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
list.of.packages = c("GSEABase","GRridge","GSA")
lapply(list.of.packages, require, character.only = TRUE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                      Load data sets                                #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
load("mirsData.RData")
attach(mirsData)

# Description:
# A deep sequencing analysis on small non-coding ribonucleic acid (miRNAseq)
# was performed on 56 samples (24 women with high-grade cervical 
# intraepithelial neoplasia (CIN3) and 32 healthy women) for the purpose of
# finding relevant screening markers for cervical cancer screening.
# The next generation sequencing analysis resulted in 2,576 transcripts. 
# The data was normalized and pre-processed, rendering 772 transcripts.
# More detail description of the data sets and the preprocessing step
# are available on the supplementary material of this following publication
# "Better diagnostic signatures from RNAseq data through use of auxiliary co-data", 
# Bioinformatics (2017).
#
# The "mirsData" object contains:
# 1. countData: a data frame of count data
# 2. transformedData: a data frame of transformed data. The countData was 
#    transformed to quasi-gausian scale.
# 3. standardizedData: a data frame of standardized data (from the transformedData).
# 4. response: a factor containing responses for control samples (n=32)
#    and samples with CIN3 (n=24)
# 5. conservation: a factor containing conservation status of each transcript.
#    The list was taken from TargetScanHuman version 7.1.
#    The conservation status of a miRNA that is divided into three classes, namely 
#    non-conserved (0), conserved across mammals (1) and broadly conserved
#    across most vertebrates (2).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                              Creating partitions                                   #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# standard deviation (sds)
Sds = apply(transformedData,1,sd)
parSds = CreatePartition(Sds,decreasing=TRUE,mingr=25,uniform=TRUE,ngroup=10)

# Conservation status
parCons = CreatePartition(conservation)

# Abundance
abundance = rowSums(countData) 
parAbund = CreatePartition(abundance,mingr=25,ngroup=10,decreasing=TRUE)

# Combining all partitions into a list
ListPartitions = list(abundance=parAbund,sds = parSds,conservation=parCons)
ListMonotone = c(TRUE,FALSE,FALSE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                    Partition selection                             #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# apply partition selection and partition ordering in order to optimize the
# performance of grridge model
# (! it may take some time)
selPar = PartitionsSelection(standardizedData, response, ListPartitions, ListMonotone, optl=NULL, innfold=NULL) 

# a list of partitions that improve the performance of GRridge model (based on cvl)
partitionsUpdate = ListPartitions[selPar$ordPar]

# a list of monotone functions for the corresponding partition
monotoneUpdate = ListMonotone[selPar$ordPar]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                 GRridge model                                      #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
grMirs = grridge(standardizedData,response,partitions=partitionsUpdate,
                 monotone=monotoneUpdate,optl=selPar$optl) 
grMirsCV = grridge.cv(grMirs,standardizedData,response)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
