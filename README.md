# Ets1 bioinformatics analysis user codes
The paper focuses on how chromatin architecture dynamics during osteoblast differentiation affect osteogenic genes, and on the mechanisms by which chromatin architecture dynamics occur in osteoblasts, and is published in AAAAA.

![Alt text](./FIG/P1.png "Ets1")

# Input files
![Alt text](./FIG/P2.png "Ets1")

* Bulk RNA-seq at DAY0, 4, 7, 10 with 4 replicates
* ChIP-seq at GM, OM about H3K27ac, H3K4me1, H3K4me3
* Cut&Run-seq at GM, OM about Ets1, p300, H3K27ac
* Bulk ATAC-seq at GM, OM with 4 replicates
* scATAC-seq about osteoblast linage cell
* Micro-C at GM, OM with 3 replicates

# Dependency

* vroom_1.6.5                              
* UpSetR_1.4.0                             
* ATACseqQC_1.26.0                         
* NucleoATACR_1.1                          
* VplotR_1.12.1                            
* TxDb.Mmusculus.UCSC.mm10.knownGene_3.10.0
* TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2  
* TxDb.Hsapiens.UCSC.hg38.knownGene_3.16.0 
* GenomicFeatures_1.54.1                   
* clusterProfiler_4.10.0                   
* ExperimentHub_2.10.0                     
* AnnotationHub_3.10.0                     
* BiocFileCache_2.10.1                     
* dbplyr_2.5.0                             
* org.Hs.eg.db_3.18.0                      
* org.Mm.eg.db_3.18.0                      
* directlabels_2024.1.21                   
* tximportData_1.30.0                      
* rtracklayer_1.62.0                       
* flashClust_1.01-2                  
