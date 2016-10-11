
DESCRIPTION:

we created a list of  platelet-expressed genes based on the joint lists of these two following studies, 

1. Gnatenko DV et al. Transcript profiling of human platelets using microarray and serial analysis of gene expression. Blood. (2003)
This study was aimed to observe gene expression transcripts in human blood platelets. Peripheral blood taken from normal healthy volunteers were hybridized its platelets in Affymetrix HG-U95Av2 (12,599 probes). Gene expression quantification (calculated as the average difference value) was applied to determine the human platelet-expressed genes. The study found 2,147 genes that were expressed
in platelets samples.  

2. Bugert P et al. Messenger RNA profiling of human platelets by microarray hybridization. ThrombHaemost. (2003)
A microarray gene expression data set containing 9,850 human genes was generated from 24 healthy Caucasians. Mean signal intensity for each and every gene was estimated, rendering 1,526 genes were found having positive hybridization signals. We included those positive genes into our platelet expressed genes list. The full list of human platelets genes is in "plt_array_Bugert P et al.xls" 

In total, the list of platelet-expressed genes contains of 1,529 genes.The genes in the RNAseq data were grouped to either “platelet-expressed genes” or “non-overlapped genes”. As it implies by the name, “platelet-expressed genes” is a group of genes in the RNAseq data that were found to be express in platelets’ samples from the two aforementioned studies. Meanwhile, the other genes that are not found in the list are grouped to “non-overlapped genes” group. 


SOURCES:
1. Gnatenko DV et al.(2003)
   http://www.bloodjournal.org/content/101/6/2285.long?sso-checked=true
2. Bugert P et al.(2003)
   https://th.schattauer.de/en/contents/archive/issue/797/manuscript/3552.html


CO-DATA PREPROCESSING: 
- Due to the data availability, we only considered the top-50 platelet-expressed genes listed in Table 2 in Gnatenko, DV et al (2003)
- We took the top 1,526 genes from "plt_array_Bugert P et al.xls" (based on the meanDens).
- The combined platelets-expressed genes from those two publications are stored at "PlateletExpressedGenesList.txt"
