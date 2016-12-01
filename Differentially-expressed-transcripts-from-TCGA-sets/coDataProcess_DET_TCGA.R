# Install the required-package from Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
library(edgeR)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 1. Extract RNAseq data sets from the TCGA portal
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#RNA sequencing data sets with these following criterions were downloaded
#(a) samples with "normal match tumor"
#(b) RNAseqV2
#see https://wiki.nci.nih.gov/display/TCGA/RNASeq+Version+2 for further details 
#(c) Data level 3 
#see https://wiki.nci.nih.gov/display/TCGA/Data+level for further details
#(d) Normalized gene files (i.e. files whose name ends with ".genes")
#see https://wiki.nci.nih.gov/display/TCGA/RNASeq+Data+Format+Specification for 
#further details


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 2. Read TCGA data in R 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# All related files for this following analysis are available at
# https://github.com/markvdwiel/GRridgeCodata/tree/master/Differentially-expressed-
# transcripts-from-TCGA-sets

# The main data files (for CRC case) can be downloaded at
# https://github.com/markvdwiel/GRridgeCodata/blob/master/Differentially-expressed-transcripts-from-TCGA-sets/CRC_UNC__IlluminaHiSeq_RNASeqV2_level3_genes.zip
# Suppose the downloaded files are stored in this following directory
setwd("~/CRC_UNC__IlluminaHiSeq_RNASeqV2_level3_genes")
listFiles = list.files()

# Read "file_manifest.txt"
# It contains a list of TCGA samples and its characteristics
# e.g. platform type, data center, array platform, level, samples' name, 
# barcode and file name
filemanifest = read.table("CRC_file_manifest.txt",header=TRUE)

# Start reading the TCGA files
a = read.table(listFiles[1],header=TRUE) #just for initialization
normalDat =  matrix(NA,dim(a),1)
cancerDat =  matrix(NA,dim(a),1)
sampleID_normal = sampleID_cancer = c()
for(i in 1:length(listFiles)){
  file = listFiles[i]
  idTemp = which(file==filemanifest[,7])
  temp0 = as.character(filemanifest[idTemp,5])
  check = unlist(strsplit(temp0,"-")) 
  
  # a file whose name has "11" ("01") as its last two digits number
  # means that such file belong to normal (cancer) tissue sample
  # see "https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode" for further details
  if(check[4]=="11"){
    temp = read.table(file,header=TRUE)
    normalDat = cbind(normalDat,temp[,2])
    sampleID_normal = c(sampleID_normal, unlist(strsplit(temp0,"-11")))
  }else if(check[4]=="01"){
    temp = read.table(file,header=TRUE)
    cancerDat = cbind(cancerDat,temp[,2])
    sampleID_cancer = c(sampleID_cancer, unlist(strsplit(temp0,"-01")))
  }
}

# Match normal and cancer samples
# Some samples do not have both normal and cancer samples)
# Such samples would be excluded
is = intersect(sampleID_normal,sampleID_cancer)
icancer = match(is, sampleID_cancer)
inormal = match(is, sampleID_normal)

normalDat = normalDat[,-1]; cancerDat = cancerDat[,-1]
normalCRCDat = normalDat[,inormal]
cancerCRCDat = cancerDat[,icancer]

# We now have two matrices from normal and cancer samples (in the same order)
# column represents individual samples
# row represents features (gene_symbol|gene_ID) 
colnames(normalCRCDat) = sampleID_normal[inormal]
colnames(cancerCRCDat) = sampleID_cancer[icancer]
rownames(normalCRCDat) = rownames(cancerCRCDat) = temp[,1]  


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 3. Data normalization
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# A function for normalization
# Normalization method: weighted trimmed mean of M-values (TMM method)
normalisation = function(dat){
  libsize <- colSums(dat)
  DGE1=DGEList(dat) 
  
  y1 <- calcNormFactors(DGE1)
  
  # Filter out genes whose rowcount is less than 5
  keep <- which(rowSums(DGE1$counts>0) >= 5) 
  
  # Relative library size; relative w.r.t geometric mean
  rellibsize <- libsize/exp(mean(log(libsize)))
  
  # Mutliplication factor per sample = normfactor*rellibsize
  nf = y1$samples[,3]*rellibsize
  
  normdat = round(sweep(dat, 2, nf, "/"))
  return(normdat)
}

# Normalize each of "normal" and "cancer" data 
normalCRCDat = normalisation(normalCRCDat)
cancerCRCDat = normalisation(cancerCRCDat)

# NOTE: We only show the TCGA data preprocessing for CRC. 
# Preprocessing for other TCGA data sets may follow the same procedure.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 4. Differentially expressed transcripts analysis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Combine normal and tumour samples
CRCDat = cbind(cancerCRCDat,normalCRCDat)
respCRC = c(replicate(ncol(cancerCRCDat),"CRC"),
            replicate(ncol(normalCRCDat),"normal_CRC"))

colonDat = cbind(cancerColonDat,normalColonDat)
respColon = c(replicate(ncol(cancerColonDat),"colon"),
              replicate(ncol(normalColonDat),"normal_colon"))

dat = cbind(CRCDat,colonDat)
resp = factor(c(respCRC,respColon))

# Start differentially expressed transcripts analysis
DGE2 = DGEList(dat)
y2 = calcNormFactors(DGE2)
design = model.matrix(~0+resp)
cont = makeContrasts(contrast="(respCRC-respnormal_CRC)-(respcolon-respnormal_colon)",
                     levels=design)
disp = estimateDisp(DEG2, design=design)
fit = glmFit(y2, design, dispersion=disp) 
lrt = glmLRT(fit,contrast=cont)
pval = lrt$table[,4]  #resulted p-values from the DETs analysis on TCGA data sets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
