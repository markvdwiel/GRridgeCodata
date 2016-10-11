
DESCRIPTION:

The processed gene expression sets were extracted from the ArrayExpress online repository. 
 
Breast cancer. The study included 94 women with breast cancer (57 patients with malignant tumor and 
37 patients with benign tumor), 31 healthy controls and 15 women who had breast cancer surgery. 
Gene expression from peripheral blood mononuclear cells samples was measured by Affymetrix Gene Chip 
HG U133 Plus 2.0. For our analysis, we selected 57 patients with malignant breast cancer and 31 normal 
controls. (ArrayExpressID: E-GEOD-27562)

Non-small cell lung carcinoma (NSCLC). A set of gene expression data from 162 subjects (82 NSCLC patients 
and 80 healthy controls), was profiled to identify differentially expressed genes between lung 
adenocarcinoma patients and healthy controls. Affymetrix Gene Chip HG U133A 2.0 was used to measure 22,277
probes from lung tissue and peripheral whole blood of the subjects. For our study, we included samples 
from whole blood only. (ArrayExpressID: E-GEOD-20189) 




SOURCES:

- Breast cancer: https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-27562/
  Please follow the aforementioned link to extract the experimental data set. 
  file:
  "sampleID_breast.txt"	  : a file contains information about sampleID and phenotype for each sample.
                            This file was obtained from "E-GEOD-27562.sdrf"
                            
- NSCLC: https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-20189/
  Please follow the aforementioned link to extract the experimental data set. 
  files: 
  "E-GEOD-20189.sdrf"	  : a full description of all samples	
  "sampleID_NSCLC.txt"	  : a file contains information about sampleID, response (class label) and 
                            smoking status for each sample. This file was obtained from "E-GEOD-20189.sdrf"




CO-DATA PROCESSING:

The "coDataProcessing_DEG.R" provides an R script to process the gene expression data sets 

The R script includes these following precedures:
- Mapping probes into the gene level and summarizing genes with multiple probe' names.
- Differentially expressed genes (DEGs) analysis by limma approach in each and every data set.
- Combining DEGs resulted from two data sets for a co-data.
