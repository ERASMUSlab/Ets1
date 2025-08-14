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

* vroom (>=1.6.5)                              
* UpSetR (>=1.4.0)                             
* ATACseqQC (>=1.26.0)                         
* NucleoATACR (>=1.1)                          
* VplotR (>=1.12.1)                            
* TxDb.Mmusculus.UCSC.mm10.knownGene (>=3.10.0)
* GenomicFeatures (>=1.54.1)                
* clusterProfiler (>=4.10.0)                  
* ExperimentHub (>=2.10.0)                  
* AnnotationHub (>=3.10.0)                    
* BiocFileCache (>=2.10.1)                    
* dbplyr (>=2.5.0)                    
* org.Hs.eg.db (>=3.18.0)                 
* org.Mm.eg.db (>=3.18.0)                     
* directlabels (>=2024.1.21)                   
