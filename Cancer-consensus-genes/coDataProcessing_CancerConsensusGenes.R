
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                   Read the downloaded cancer consensus list                      #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
COSMIC = read.csv("Census_allWed Nov 25 11_56_18 2015.csv")



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#  Re-group the somatic consensus genes based on the type of cancers of interest   #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# The listed genes were next used to group genes in the primary data set into three 
# groups, namely (1) cancer specific genes, i.e. genes that are known to mutate in
# cancer diseases of interest (2) cancer genes, i.e. known mutated genes in the 
# other than the cancer specific genes; and (3) genes that were not stated as a
# cancer mutated gene by the COSMIC list. 
# This section provides a preprocessing step to re-group the cancer somatic genes
# into group (1) and group (2).
# Example: Lung and Pancreas cancer case
idLung = idPancreas = idBoth = c()
# idLung/idPancreas/idBoth: a vector containing indeces of mutated genes in 
# lung cancer/pancreatic cancer/both cancers
for(i in 1:dim(COSMIC)[1]){
  temp = as.character(COSMIC[i,8])
  temp2 = unlist(strsplit(temp,", "))
  
  # As one gene possibly cause more than one cancer,
  # we check which mutated genes in the spesific cancers of interest,
  # in this case, lung and pancreas cancer
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

# Re-group the somatic cancer genes into (1) NSCLCpancreatic and (2) other cancers
idAll = unique(c(idLung,idPancreas,idBoth))
cosmic =  replicate(nrow(COSMIC),"other cancers")
cosmic[idAll] = "NSCLCpancreatic"
names(cosmic) = COSMIC[,1]
table(cosmic)

# NOTE:
# 
# cosmic
# NSCLCpancreatic   other cancers 
# 42             530
# 
# it means that there are 42 genes (out of 572 cancer somatic genes) that are found
# to be mutated both in NSCLC and in Pancreatic cancer.