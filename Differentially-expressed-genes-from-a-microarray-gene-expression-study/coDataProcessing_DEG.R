
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                   A list of required-packages for this analysis                  #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
list.of.packages <- c("affy","genefilter","AffyExpress","hgu133plus2.db","hgu133a.db")
lapply(list.of.packages, require, character.only = TRUE)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                      A FUNCTION                                  #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# A function to:
# (1) map probes to the gene level and 
# (2) summarize genes with multiple names of probe
# For a gene with multiple names of probe, we took a median of their expression values.
# Please run this function first prior to the differentially expressed genes analysis

duplicationAffy = function(version,data,n){
  ids=rownames(data)
  GENENAME = get(paste(version,"GENENAME", sep=""))
  ENTREZID = get(paste(version,"ENTREZID", sep=""))
  GENESYMBOL = get(paste(version,"SYMBOL", sep="")) 
  
  genename = GENENAME[ids]
  genesymbol = GENESYMBOL[ids]
  probe.id = toTable(genename)$probe_id 
  entrez = ENTREZID[ids]
  gene.id = toTable(entrez)$gene_id
  genesymbol.id = toTable(genesymbol)$symbol
  
  # genename summary
  gn = merge(toTable(entrez), toTable(genesymbol))
  
  # unknown probeset's names
  mp = ids[is.na(pmatch(ids,probe.id))] ; length(mp) #missing probesets
  mp.mat = data[which(is.na(pmatch(ids,probe.id))=="TRUE"),] 
  
  # Newdata: expressions on  known probesets
  data.temp = data[which(is.na(pmatch(ids,probe.id))=="FALSE"),]
  
  # unique genes
  gene.uniq = unique(gn[,3]); p2=length(gene.uniq);
  
  newdata = matrix(NA,p2,n) ; 
  for(i in 1:p2){
    idx = which(gn[,3] == gene.uniq[i])
    for(j in 1:n){newdata[i,j] = median(data.temp[idx,j])}
  }
  rownames(newdata) = gene.uniq
  
  sum = matrix(NA,3,1)
  rownames(sum)=c("unknown probesets","known probesets","unique genes")
  sum[1,1] = length(mp)
  sum[2,1] = dim(gn)[1]
  sum[3,1] = as.numeric(p2)
  print(sum)
  
  summary=list()
  summary[[1]]=newdata #unduplicated genes
  summary[[2]]=gn #known probesets 
  summary[[3]]= mp.mat #unknown probesets
  
  return(summary)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#               DIFFERENTIALLY EXPRESSED GENES ANALYSIS FOR BREAST CANCER DATA               #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# (1) A preprocessed file for each sample is available at 
#     https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-27562/samples/
#
# (2) A file containing more information for each sample 
#     is available at 
#     https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-27562/E-GEOD-27562.sdrf.txt
#     We summarized the file as such it is ready to be used for this analysis
#     It is available at
#     https://github.com/markvdwiel/GRridgeCodata/tree/master/Differentially
#     -expressed-genes-from-a-microarray-gene-expression-study
#    (i.e. "sampleID_breast.txt")

matID = read.table("sampleID_Breast.txt",header=TRUE)
# "sampleID_Breast.txt" is a file containing more information for each sample, which can be
# obtained from the aforementioned link. We took information about the sampleID and phenotype 
# only and then save such information to the new "sampleID_Breast.txt" file.
# Additionally, we only took "malignant" samples and "normal" controls.
sampleID = matID[,1]
respBreast = as.factor(matID[,2]) 
nrespBreast = length(respBreast)

matBreast = matrix(NA,54675,1) #54675 is the initial number of probes from the array
for(i in 1: nrespBreast){
  temp = as.matrix(read.table(paste(sampleID[i],"_sample_table.txt",sep=""),header=TRUE))
  matBreast = cbind(matBreast,as.numeric(temp[,2]))
}
matBreast = matBreast [,-1]
colnames(matBreast) = respBreast
rownames(matBreast) = temp[,1]
dim(matBreast)

# summarize duplications 
matBreast2 = duplicationAffy("hgu133plus2", matBreast, nrespBreast)[[1]]

design = model.matrix(~0+respBreast)
# Fit linear model
fit = lmFit(matBreast2, design);   
cont = makeContrasts(respBreastMalignant-respBreastNormal,levels=design)
cont.fit = contrasts.fit(fit,cont)
fite = eBayes(cont.fit)

resBreast = topTable(fite, number=nrow(matBreast2), adjust="fdr")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                   DIFFERENTIALLY EXPRESSED GENES ANALYSIS FOR NSCLC DATA                   #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#(1) A preprocessed file for each sample is available at 
#    https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-20189/samples/
# 
#(2) A file containing more information for each sample ("sampleID_NSCLC.txt")
#    is available at 
#    https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-20189/E-GEOD-20189.sdrf.txt
#    We summarized the file as such it is ready to be used for this analysis
#    It is available at
#    https://github.com/markvdwiel/GRridgeCodata/tree/master/Differentially- 
#    expressed-genes-from-a-microarray-gene-expression-study
#    (i.e. "sampleID_NSCLC.txt")

matID2 = read.table("sampleID_NSCLC.txt",header=TRUE)
# "sampleID_NSCLC.txt" is a file containing more information for each sample, which can be
# obtained from the aforementioned link. We took information about the sampleID and phenotype 
# only and the save such information to the new "sampleID_NSCLC.txt" file.
# We took information about the sampleID, phenotype and smoking status.
# and then saved such information to the new "sampleID_NSCLC.txt" file.

sampleID2 = matID2[,1]
respNSCLC = as.factor(matID2[,2]) #case:NSCLC patients; control:HC
nrespNSCLC = length(respNSCLC)
smooke = as.factor(matID2[,3])

matNSCLC = matrix(NA,22277,1) #22277 is the initial number of probes from the array
for(i in 1: nrespNSCLC){
  temp = as.matrix(read.table(paste(sampleID2[i],"_sample_table.txt",sep=""),header=TRUE))
  matNSCLC = cbind(matNSCLC,as.numeric(temp[,2]))
}
matNSCLC = matNSCLC[,-1]
colnames(matNSCLC) = respNSCLC
rownames(matNSCLC) = temp[,1]

# summarize duplications
matNSCLC2 = duplicationAffy("hgu133a",matNSCLC,nrespNSCLC)[[1]]

# adjust limma model by "smoking" status
design = model.matrix(~0+respNSCLC+smooke)
fit = lmFit(matNSCLC2, design);   # Fit linear matrix
cont = makeContrasts(respNSCLCCase-respNSCLCControl,levels=design)
cont.fit = contrasts.fit(fit,cont)
fite = eBayes(cont.fit)

respNSCLC = topTable(fite, number=nrow(matNSCLC2), adjust="fdr")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                            Combine the two DEG analysis results                            #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
FDR = 0.05
idDEGNSCLC = which(respNSCLC[,5]<FDR)
degNSCLC=rownames(respNSCLC)[idDEGNSCLC]
idDEGBreast = which(resBreast[,5]<FDR)
degBreast=rownames(resBreast)[idDEGBreast]
degBoth = intersect(degBreast,degNSCLC)

# Assign the DEG information to the primary data set
# "genesWurdinger" is an object (vector), containing the gene symbol from the mRNAseq data set
# help(dataWurdinger) shows how to get the "genesWurdinger" object. 
exprSet = matrix("noOverlap",1,length(genesWurdinger))
# DEGs from the Breast cancer data set
im = intersect(genesWurdinger,degBreast)
is = match(im, genesWurdinger)
exprSet[is] = "deg_breast"
# DEGs from the NSCLC data set
im = intersect(genesWurdinger,degNSCLC)
is = match(im, genesWurdinger)
exprSet[is]="deg_NSCLC"
# Genes that are stated as differentially expressed genes by both data sets
im = intersect(genesWurdinger,degBoth)
is = match(im, genesWurdinger)
exprSet[is]="deg_both"
exprSet = as.factor(exprSet)
table(exprSet)


# exprSet is an object that can be used to create a partition for the purpose of GRridge modeling
# see help("CreatePartitions") for more details on how to create a partition in the GRridge package
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
