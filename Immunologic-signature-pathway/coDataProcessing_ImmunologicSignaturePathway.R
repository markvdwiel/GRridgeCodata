# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                    Get information from the GSA databases                        #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# IMPORTANT NOTE: This part is not incorporated in the GRridge model
# It is aimed to observe the details of the downloaded object from the gene set
# enrichment analysis (GSEA). See help(GSA.read.gmt) for further detail.
library(GSA)  
gsaImmun.sym = GSA.read.gmt("c7.all.v5.1.symbols.gmt")
genesetsImmun.sym = gsaImmun.sym[[1]]
genesetNamesImmun.sym = gsaImmun.sym[[2]]
genesetDescImmun.sym = gsaImmun.sym[[3]]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                              CO-DATA PROCESSING                                  #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Exporting a GMT file to be used further in a GRridge model 
library(GSEABase)
ImmunSym = getGmt("c7.all.v5.1.symbols.gmt")

# Creating partitions based on pathways information (e.g. GSEA object)
# Some variables may belong to more than one groups (gene sets)
# The argument minlen=25 implies the minimum number of members of a gene set
# If remain=TRUE, gene sets with less than 25 members are grouped to the
# "remainder" group
# "genesWurdinger" is an object containing genes in the data set in the mRNA example
# see help(matchGeneSets) for detail information of arguments in the function.
gseImmun = matchGeneSets(genesWurdinger,ImmunSym,minlen=25,remain=TRUE)

# Regrouping gene sets by hierarchical clustering analysis
# The number of gene sets from the GSEA database is relatively too high to be used 
# in the GRridge model. Here, the initial gene sets are re-grouped into maxGroups=5,
# using information from the primary data set.
# "datStdWurdinger" is an object containing standardized RNAseq data
gseImmun_newGroups = mergeGroups(highdimdata= datStdWurdinger, initGroups =gseImmun, maxGroups=5);

# Extracting indeces of new groups
gseImmun2 = gseImmun_newGroupst$newGroups

# Members of the new groups
newClustMembers = gseImmun_newGroups$newGroupMembers
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

