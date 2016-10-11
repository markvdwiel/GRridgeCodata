DESCRIPTION:
The cancer genome atlas (TCGA) is a public database that contains a considerable number of samples
used in cancer experiments and it covers many types of omics data sets (e.g. DNA methylation, 
copy number, array-based expression and mRNA sequencing) and clinical information. 
RNASeq data sets from breast cancer (n=110), colon adenocarcinoma (n=26) and lung adenocarcinoma (n=58)
were available at the public repository (last checked on June 3, 2016). 
We aim to extract p-values from those data sets to be used further as a co-data in GRridge modeling, 
in a binary classification case that include those cancer types. 
Binary classification cases included those cancer types may use p-values information as a co-data. 

SOURCES:
TCGA website.
RNA sequencing data sets with these following criterion were downloaded:
(a) samples with "normal match tumor"
(b) RNAseqV2 (see https://wiki.nci.nih.gov/display/TCGA/RNASeq+Version+2)
(c) Data level 3 (see https://wiki.nci.nih.gov/display/TCGA/Data+level)
(d) Normalized gene files (i.e. files whose name ends with “.genes”)
    (see https://wiki.nci.nih.gov/display/TCGA/RNASeq+Data+Format+Specification)

files:
- Breast cancer
  "BreastCancer_UNC__IlluminaHiSeq_RNASeqV2_level3_genes": contains normalized file for each sample
  "BreastCancer_file_manifest.txt": a file contains detail description of each sample
- Colorectal cancer (CRC)
  "CRC_UNC__IlluminaHiSeq_RNASeqV2_level3_genes": contains normalized file for each sample  
  "CRC_file_manifest.txt": a file contains detail description of each sample
- Non-small cell lung carcinoma (NSCLC)
  "NSCLC_UNC__IlluminaHiSeq_RNASeqV2_level3_genes": contains normalized file for each sample  
  "NSCLC_file_manifest.txt": a file contains detail description of each sample


CO-DATA PROCESSING:
The "coDataProcess_DET_TCGA.R" provides an R script to process the TCGA data sets as a co-data for GRridge
modeling in the RNASeq data set.
 
The R script includes these following precedures:
- Read the downloaded TCGA data sets in the R environment
- Normalize the data sets
- Perform differentially expressed transcripts (P values from this analysis may be used as a co-data)
