# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#          1. Extract information of the cancer consensus genes               #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# source: http://cancer.sanger.ac.uk/census
# The cancer consensus genes used in this study was stamped on November 25, 2015 
# It contained of 572 cancer somatic genes. 
# The file is also available at:
# https://github.com/markvdwiel/GRridgeCodata/tree/master/Cancer-consensus-genes
# "Census_allWed Nov 25 11_56_18 2015.csv"

datCosmic = read.csv("Census_allWed Nov 25 11_56_18 2015.csv")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#               2. Re-group the somatic consensus genes                       #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# The grouping is based on the type of cancers of interest
idLung = idPancreas = idBoth = c()
# idLung/idPancreas/idBoth: a vector containing indices of mutated genes in 
# lung cancer/pancreatic cancer/both cancers
for(i in 1:dim(datCosmic)[1]){
  temp = as.character(datCosmic[i,8])
  temp2 = unlist(strsplit(temp,", "))
  
  # As one gene possibly causes more than one cancer,
  # we check which mutated genes in the spesific cancers of interest 
  #(i.e. NSCLC and pancreatic cancer)
  checkA1 = "lung" %in% temp2
  checkA2 = "NSCLC" %in% temp2
  checkB1 = "pancreas" %in% temp2
  checkB2 = "pancreatic" %in% temp2
  
  if(checkA1 ==TRUE || checkA2==TRUE){idLung=c(idLung,i)}
  if(checkB1 ==TRUE || checkB2==TRUE) {idPancreas=c(idPancreas,i)}
  if((checkA1 ==TRUE && checkB1==TRUE) || (checkA1 ==TRUE && checkB2==TRUE) ||
     (checkA2 ==TRUE && checkB1==TRUE) || (checkA2 ==TRUE && checkB2==TRUE)) {idBoth=c(idBoth,i)}
  else{next}
}

# Re-group the somatic cancer types into (1) NSCLCpancreatic and (2) other cancer
idAll = unique(c(idLung,idPancreas,idBoth))
cosmic =  replicate(nrow(datCosmic),"other cancers")
cosmic[idAll] = "NSCLCpancreatic"
table(cosmic)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#       3. Assign the somatic cancer genes information to each of feature     #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# "genesWurdinger" is an object (vector), 
# containing the gene symbol from the mRNAseq data set
# help(dataWurdinger) shows how to get the "genesWurdinger" object.
is = intersect(datCosmic[,1], genesWurdinger)
ig = match(is, genesWurdinger)
ic = match(is,datCosmic[,1])
cosmic2 = replicate(length(genesWurdinger),"noncosmic")
cosmic2[ig] = cosmic[ic]
cosmic2 = as.factor(cosmic2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# NOTE:
# 
# cosmic
# NSCLCpancreatic   other cancers 
# 42             530
# 
# it means that there are 42 genes (out of 572 cancer somatic genes) that are found
# to be mutated both in NSCLC and in Pancreatic cancer.