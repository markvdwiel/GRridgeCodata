
# In this R script, we give an example for the co-data preparation in the case 
# of breast cancer and colon adenocarcinoma in the RNAseq data
# Preprocessing for other TCGA data sets may follow the same procedure.


library(edgeR)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                           Read TCGA data in R                               # 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Here, we only show how to read TCGA data in R in the Breast Cancer case

setwd("Breastcancer_UNC__IlluminaHiSeq_RNASeqV2/Level_3")
listFiles = list.files()

# Read "file_manifest.txt"
# The file contains a list of TCGA samples and its characteristics
# e.g. platform type, data center, array platform, level, samples' name, 
# barcode and file name
filemanifest = read.table("BreastCancer_file_manifest.txt",header=TRUE)

# Start reading the TCGA files
a = read.table(listFiles[1],header=TRUE) #for initialization
normalDat =  matrix(NA,dim(a),1)
cancerDat =  matrix(NA,dim(a),1)
sampleID_normal = sampleID_cancer = c()
for(i in 1:length(listFiles)){
  file = listFiles[i]
  idTemp = which(file==filemanifest[,7])
  temp0 = as.character(filemanifest[idTemp,5])
  check = unlist(strsplit(temp0,"-")) 
  
  # A file whose name has "11" ("01") as its last two digits number
  # means that such file belong to normal (cancer) tissue sample
  # see "https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode" for further detail
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
# note: some samples do not have both normal and cancer samples
# Such samples would be excluded
is = intersect(sampleID_normal,sampleID_cancer)
icancer = match(is, sampleID_cancer)
inormal = match(is, sampleID_normal)

normalDat = normalDat[,-1]; cancerDat = cancerDat[,-1]
normalBreastDat = normalDat[,inormal]
cancerBreastDat = cancerDat[,icancer]

# We now have two matrices from normal and cancer samples (in the same order)
# column represents individual samples
# row represents features (gene_symbol|gene_ID) 
colnames(normalBreastDat) = sampleID_normal[inorm]
colnames(cancerBreastDat) = sampleID_cancer[ic]
rownames(normalBreastDat) = rownames(cancerBreastDat) = temp[,1]  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                  Normalization                              #
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
# Apply the "normalization" function
normalBreastDat = normalisation(normalBreastDat)
cancerBreastDat = normalisation(cancerBreastDat)

# NOTE: the same normalization procedure for other TCGA sets can follow the
# the same procedure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#             Differentially expressed transcripts analysis                   #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Combine normal and tumour samples in each cancer case
BreastDat = cbind(cancerBreastDat,normalBreastDat)
respNsclc = c(replicate(ncol(cancerBreastDat),"nsclc"),
              replicate(ncol(normalBreastDat),"normal_nsclc"))

colonDat = cbind(cancerColonDat,normalCOlonDat)
respColon = c(replicate(ncol(cancerColonDat),"colon"),
              replicate(ncol(normalCOlonDat),"normal_colon"))

# Create a new matrix for combined data sets
dat = cbind(nsclcDat,colonDat)
resp = factor(c(respNsclc,respColon))

# Start differentially expressed transcripts analysis
DGE2 = DGEList(dat)
y2 = calcNormFactors(DGE2)
design = model.matrix(~0+resp)
cont = makeContrasts(contrast="(respnsclc-respnormal_nsclc)-(respcolon-respnormal_colon)",levels=design)
disp = estimateDisp(DEG2, design=desigh)
fit = glmFit(y2, design, dispersion=disp) 
lrt = glmLRT(fit,contrast=cont)
pval = lrt$table[,4]  #resulted p-values from the DETs analysis on TCGA data sets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #