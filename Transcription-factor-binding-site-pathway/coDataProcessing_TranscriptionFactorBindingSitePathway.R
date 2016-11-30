
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                    Get information from the GSA databases                        #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# IMPORTANT NOTE: This part is not incorporated in the GRridge model
# It is aimed to observe the details of the downloaded object from the gene set
# enrichment analysis (GSEA). See help(GSA.read.gmt) for further detail.
library(GSA)  
gsaTF.sym = GSA.read.gmt("c3.tft.v5.0.symbols.gmt")
genesetsTF.sym = gsaTF.sym[[1]]
genesetNamesTF.sym = gsaTF.sym[[2]]
genesetDescTF.sym = gsaTF.sym[[3]]



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                              CO-DATA PROCESSING                                  #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Exporting a GMT file to be used further in a GRridge model 
library(GSEABase)
TFsym = getGmt("c3.tft.v5.0.symbols.gmt")

# Creating partitions based on pathways information (e.g. GSEA object)
# Some variables may belong to more than one groups (gene sets).
# The argument minlen=25 implies the minimum number of members in a gene set
# If remain=TRUE, gene sets with less than 25 members are grouped to the
# "remainder" group.
# help(dataWurdinger): shows on how to get the "genesWurdinger" object from the GRridge package
# help(matchGeneSets): for detail information of the "matchGeneSets" function.
gseTF = matchGeneSets(genesWurdinger,TFsym,minlen=25,remain=TRUE)

# Regrouping gene sets by hierarchical clustering analysis
# The number of gene sets from the GSEA database is relatively too high to be used 
# in the GRridge model. Here, the initial gene sets are re-grouped into maxGroups=5, 
# using information from the primary data set.
gseTF_newClust = mergeGroups(highdimdata= datStd, initGroups =gseTF, maxGroups=5);

# Extracting indeces of new groups
gseTF2 = gseTF_newClust$newGroups

# Members of the new groups
newClustMembers = gseTF_newClust$newGroupMembers
