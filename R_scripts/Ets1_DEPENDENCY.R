#install.packages("ggplot2") 
#install.packages("dplyr")
#install.packages("uwot") 
#install.packages("UpSetR") 
#install.packages("vroom") 
#install.packages("pheatmap") 
#install.packages("forcats") 
#install.packages("gridExtra") 
#install.packages("RColorBrewer")

#remotes::install_github("GreenleafLab/NucleoATACR")
#remotes::install_github("robinweide/GENOVA")

#micromamba install -c bioconda -c conda-forge -c r r-data.table bioconductor-vplotr bioconductor-clusterprofiler bioconductor-deseq2

#BiocManager::install("org.Mm.eg.db")
#BiocManager::install("org.Hs.eg.db")


options(scipen=999)
set.seed(7777)

library(uwot)
library(DOSE)
library(vroom)
library(dplyr)
library(DESeq2)
library(UpSetR)
library(VplotR)
library(GENOVA)
library(ggplot2)
library(forcats)
library(pheatmap)
library(gridExtra)
library(data.table)
library(NucleoATACR)
library(RColorBrewer)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)

