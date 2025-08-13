[scLIKE_ver2.md](https://github.com/user-attachments/files/21748227/scLIKE_ver2.md)
```R
options(scipen=999)

# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("DESeq2")
#install.packages("DESeq2")

library(DT)
library(DOSE)
library(vroom)
library(WGCNA)
library(readr)
library(dplyr)
library(tidyr)
library(qqman)
library(GO.db)
library(caret)
library(fgsea)
library(VplotR)
library(tibble)
library(ChIPQC)
library(GOplot)
library(DESeq2)
library(UpSetR)
library(viridis)
library(ggrepel)
library(forcats)
library(ggplot2)
library(biomaRt)
library(msigdbr)
library(viridis)
library(CEMiTool)
library(pheatmap)
library(tximport)
library(corrplot)
library(GOSemSim)
library(circlize)
library(ATACseqQC)
library(patchwork)
library(tidyverse)
library(DMRcaller)
library(networkD3)
library(hrbrthemes)
library(ChIPseeker)
library(ReactomePA)
library(enrichplot)
library(data.table)
library(multiWGCNA)
library(flashClust)
library(NucleoATACR)
library(rtracklayer)
library(tximportData)
library(directlabels)
library(BiocParallel)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(RColorBrewer)
library(ExperimentHub)
library(GenomicRanges)
library(dynamicTreeCut)
library(clusterProfiler)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
```


```R
UMAP_viewer = function(inputPATH){

UMAP = read.table(inputPATH,sep = "\t", header=T)

options(repr.plot.width = 8, repr.plot.height = 8, repr.plot.res = 200, repr.plot.pointsize = 30)

color_palette <- c("darkblue", "#B71C1C","#283593", "#3949AB", "#5C6BC0", "#9FA8DA",
                   "#757575", "#FFCDD2", "#E57373", "#E53935")

ggplot(UMAP, aes(X, Y, color = factor(DIFF_Q))) +
  geom_point(size = 1, alpha = 1) +
  scale_color_manual(values = color_palette) +
  theme_classic(base_size = 10) +
  ggtitle(NULL) +
  xlab("UMAP1") +
  ylab("UMAP2") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 25, color = "darkblue"),
    axis.text.x = element_text(face = "bold", size = 20),
    axis.text.y = element_text(face = "bold", size = 20),
    axis.title.x = element_text(face = "bold", size = 20, color = "darkblue"),
    axis.title.y = element_text(face = "bold", size = 20, color = "darkblue"),
    legend.position = "none",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(face = "bold", size = 7)
  )
    }
```
mm10_chr_size = read.table("/data3/psg/NGS_2025/4DN/HiC_data/mm10.chrom.sizes",sep = "\t", header=T)
mm10_chr_size$bin = ceiling(mm10_chr_size[,2]/100000)

chrlist = mm10_chr_size[,1]
chrlist = as.data.frame(chrlist)
colnames(chrlist) = "chr"

for(i in 1:19){
chrlist_DF = mm10_chr_size[which(mm10_chr_size$chr == chrlist[i,1]),]
options(scipen=999)
scale = 100000
scale = as.integer(scale)
end = chrlist_DF[1,3]

s = 1
bin_name = paste(chrlist[i,1],(((s-1)*scale)+1),s*scale)
bin_name = as.data.frame(bin_name)
colnames(bin_name) = "bin"

for(s in 2:end){
bin_name_frag = paste(chrlist[i,1],(((s-1)*scale)+1),s*scale)
bin_name_frag = as.data.frame(bin_name_frag)
colnames(bin_name_frag) = "bin"
bin_name = rbind(bin_name,bin_name_frag)
    }

m = 1
bin_DF = cbind(bin_name,bin_name[m,1])
colnames(bin_DF) = c("X","Y")

for(m in 2:nrow(bin_name)){
bin_DF_frag = cbind(bin_name,bin_name[m,1])
colnames(bin_DF_frag) = c("X","Y")
bin_DF = rbind(bin_DF,bin_DF_frag)
    }

path = paste("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr",i,"_binDF.txt",sep="")
write.table(bin_DF, file = path, col.names=F, row.names=F, quote=F, sep="\t")
    }
chr1_binDF = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr1_binDF.txt",sep = "\t", header=F)
chr2_binDF = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr2_binDF.txt",sep = "\t", header=F)
chr3_binDF = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr3_binDF.txt",sep = "\t", header=F)
chr4_binDF = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr4_binDF.txt",sep = "\t", header=F)
chr5_binDF = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr5_binDF.txt",sep = "\t", header=F)
chr6_binDF = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr6_binDF.txt",sep = "\t", header=F)
chr7_binDF = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr7_binDF.txt",sep = "\t", header=F)
chr8_binDF = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr8_binDF.txt",sep = "\t", header=F)
chr9_binDF = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr9_binDF.txt",sep = "\t", header=F)
chr10_binDF = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr10_binDF.txt",sep = "\t", header=F)
chr11_binDF = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr11_binDF.txt",sep = "\t", header=F)
chr12_binDF = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr12_binDF.txt",sep = "\t", header=F)
chr13_binDF = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr13_binDF.txt",sep = "\t", header=F)
chr14_binDF = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr14_binDF.txt",sep = "\t", header=F)
chr15_binDF = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr15_binDF.txt",sep = "\t", header=F)
chr16_binDF = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr16_binDF.txt",sep = "\t", header=F)
chr17_binDF = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr17_binDF.txt",sep = "\t", header=F)
chr18_binDF = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr18_binDF.txt",sep = "\t", header=F)
chr19_binDF = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr19_binDF.txt",sep = "\t", header=F)

binDF_col = rbind(chr1_binDF,
                    chr2_binDF,
                    chr3_binDF,
                    chr4_binDF,
                    chr5_binDF,
                    chr6_binDF,
                    chr7_binDF,
                    chr8_binDF,
                    chr9_binDF,
                    chr10_binDF,
                    chr11_binDF,
                    chr12_binDF,
                    chr13_binDF,
                    chr14_binDF,
                    chr15_binDF,
                    chr16_binDF,
                    chr17_binDF,
                    chr18_binDF,
                    chr19_binDF)

write.table(binDF_col, file = "/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/FULL_binDF.txt", col.names=F, row.names=F, quote=F, sep="\t")mm10_chr_size = read.table("/data3/psg/NGS_2025/4DN/HiC_data/mm10.chrom.sizes",sep = "\t", header=T)
mm10_chr_size$bin = ceiling(mm10_chr_size[,2]/100000)

chrlist = mm10_chr_size[,1]
chrlist = as.data.frame(chrlist)
colnames(chrlist) = "chr"

i = 1
for(i in 1:nrow(chrlist)){
chrlist_DF = mm10_chr_size[which(mm10_chr_size$chr == chrlist[i,1]),]
end = chrlist_DF[1,3]

chrlist_DF = paste(chrlist[i,1],1,sep="_")

for(s in 2:end){
    chrlist_DF_frag = paste(chrlist[i,1],s,sep="_")
    chrlist_DF = c(chrlist_DF,chrlist_DF_frag)
    }

chrlist_DF = as.data.frame(chrlist_DF)

path = paste("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr",i,"_binNAME.txt",sep="")
write.table(chrlist_DF, file = path, col.names=F, row.names=F, quote=F, sep="\t")
    }


chr1_binNAME = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr1_binNAME.txt",sep = "\t", header=F)
chr2_binNAME = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr2_binNAME.txt",sep = "\t", header=F)
chr3_binNAME = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr3_binNAME.txt",sep = "\t", header=F)
chr4_binNAME = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr4_binNAME.txt",sep = "\t", header=F)
chr5_binNAME = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr5_binNAME.txt",sep = "\t", header=F)
chr6_binNAME = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr6_binNAME.txt",sep = "\t", header=F)
chr7_binNAME = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr7_binNAME.txt",sep = "\t", header=F)
chr8_binNAME = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr8_binNAME.txt",sep = "\t", header=F)
chr9_binNAME = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr9_binNAME.txt",sep = "\t", header=F)
chr10_binNAME = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr10_binNAME.txt",sep = "\t", header=F)
chr11_binNAME = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr11_binNAME.txt",sep = "\t", header=F)
chr12_binNAME = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr12_binNAME.txt",sep = "\t", header=F)
chr13_binNAME = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr13_binNAME.txt",sep = "\t", header=F)
chr14_binNAME = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr14_binNAME.txt",sep = "\t", header=F)
chr15_binNAME = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr15_binNAME.txt",sep = "\t", header=F)
chr16_binNAME = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr16_binNAME.txt",sep = "\t", header=F)
chr17_binNAME = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr17_binNAME.txt",sep = "\t", header=F)
chr18_binNAME = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr18_binNAME.txt",sep = "\t", header=F)
chr19_binNAME = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/chr19_binNAME.txt",sep = "\t", header=F)

binNAME_col = rbind(chr1_binNAME,
                    chr2_binNAME,
                    chr3_binNAME,
                    chr4_binNAME,
                    chr5_binNAME,
                    chr6_binNAME,
                    chr7_binNAME,
                    chr8_binNAME,
                    chr9_binNAME,
                    chr10_binNAME,
                    chr11_binNAME,
                    chr12_binNAME,
                    chr13_binNAME,
                    chr14_binNAME,
                    chr15_binNAME,
                    chr16_binNAME,
                    chr17_binNAME,
                    chr18_binNAME,
                    chr19_binNAME)

write.table(binNAME_col, file = "/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/FULL_binNAME.txt", col.names=F, row.names=F, quote=F, sep="\t")FULL_binNAME = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/FULL_binNAME.txt",sep = "\t", header=F)
FULL_binDF = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/FULL_binDF.txt",sep = "\t", header=F)

MAT_col = FULL_binNAME[,1]
MAT_row = FULL_binDF[,1]
MAT_row = as.data.frame(MAT_row)
MAT_row = unique(MAT_row)

mat = matrix(0, nrow = nrow(FULL_binNAME), ncol = nrow(FULL_binNAME))
MATRIX = as.data.frame(mat)

colnames(MATRIX) = MAT_col
rownames(MATRIX) = MAT_row[,1]FULL_binDF = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DF/FULL_binDF.txt",sep = "\t", header=F)
FULL_binDF = FULL_binDF[,1]
FULL_binDF = as.data.frame(FULL_binDF)
colnames(FULL_binDF) = c("row")

FULL_binDF = unique(FULL_binDF)

FULL_MATRIX_DF = cbind(FULL_binDF,FULL_binDF[1,1])
FULL_MATRIX_DF = as.data.frame(FULL_MATRIX_DF)
colnames(FULL_MATRIX_DF) = c("row","col")

for(i in 2:nrow(FULL_binDF)){
    FULL_MATRIX_DF_frag = cbind(FULL_binDF,FULL_binDF[i,1])
    FULL_MATRIX_DF_frag = as.data.frame(FULL_MATRIX_DF_frag)
    colnames(FULL_MATRIX_DF_frag) = c("row","col")
    FULL_MATRIX_DF = rbind(FULL_MATRIX_DF,FULL_MATRIX_DF_frag)
    }

write.table(FULL_MATRIX_DF, file = "/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/FULL_MATRIX_DF.txt", col.names=T, row.names=F, quote=F, sep="\t")args <- commandArgs(trailingOnly = TRUE)
para <- as.integer(args[1])

library(dplyr)
options(scipen=999)

MATRIX_maker = function(para){

str = para -999

MATRIX_backbone = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_backbone.txt",sep = "\t", header=T)

DAY0_WT_100kb = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/cool/DAY0_WT_100kb.txt",sep = "\t", header=F)
DAY0_WT_100kb$FULLlabel = paste(DAY0_WT_100kb[,1],DAY0_WT_100kb[,2],DAY0_WT_100kb[,3],
                                DAY0_WT_100kb[,4],DAY0_WT_100kb[,5],DAY0_WT_100kb[,6])

DAY0_WT_100kb = cbind(DAY0_WT_100kb[,8],DAY0_WT_100kb[,7])
colnames(DAY0_WT_100kb) = c("FULLlabel","interaction")
DAY0_WT_100kb = as.data.frame(DAY0_WT_100kb)

FULL_bin_region = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/FULL_bin_region.txt",sep = "\t", header=F)

FULL_bin_region = FULL_bin_region[!grepl("random",FULL_bin_region[,1]),]
FULL_bin_region = FULL_bin_region[!grepl("chrUn",FULL_bin_region[,1]),]
FULL_bin_region = FULL_bin_region[!grepl("chrX",FULL_bin_region[,1]),]
FULL_bin_region = FULL_bin_region[!grepl("chrY",FULL_bin_region[,1]),]
FULL_bin_region = FULL_bin_region[!grepl("chrM",FULL_bin_region[,1]),]

FULL_bin_region[,2] = FULL_bin_region[,2] -1
FULL_bin_region$label1 = paste(FULL_bin_region[,1],FULL_bin_region[,2],FULL_bin_region[,3])
FULL_bin_region$label2 = paste(FULL_bin_region[,1],FULL_bin_region[,2],FULL_bin_region[,3])

#for(i in 1:nrow(FULL_bin_region)){
for(i in str:para){

frag = as.data.frame(cbind(FULL_bin_region[i,4],FULL_bin_region[,5]))
frag$FULLlabel = paste(frag[,1],frag[,2])

frag = frag[,3]
frag = as.data.frame(frag)
colnames(frag) = "FULLlabel"

frag = left_join(frag,DAY0_WT_100kb,by="FULLlabel")
frag[is.na(frag)] = 0

MATRIX_backbone[i,] = frag[,2]
    }

path = paste( "/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/","MATRIX_DAY0","_",para,".txt",sep="")

write.table(MATRIX_backbone,
            file = path,
            col.names=T, row.names=T, quote=F, sep="\t")
    }

MATRIX_maker(para)MATRIX_DAY0_1000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_1000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY0_2000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_2000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY0_3000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_3000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY0_4000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_4000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY0_5000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_5000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY0_6000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_6000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY0_7000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_7000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY0_8000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_8000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY0_9000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_9000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY0_10000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_10000_header.txt", delim = "\t",show_col_types = FALSE))

rowname = MATRIX_DAY0_1000[,1]

MATRIX_DAY0_1000_fit = MATRIX_DAY0_1000[,-1]
MATRIX_DAY0_2000_fit = MATRIX_DAY0_2000[,-1]
MATRIX_DAY0_3000_fit = MATRIX_DAY0_3000[,-1]
MATRIX_DAY0_4000_fit = MATRIX_DAY0_4000[,-1]
MATRIX_DAY0_5000_fit = MATRIX_DAY0_5000[,-1]
MATRIX_DAY0_6000_fit = MATRIX_DAY0_6000[,-1]
MATRIX_DAY0_7000_fit = MATRIX_DAY0_7000[,-1]
MATRIX_DAY0_8000_fit = MATRIX_DAY0_8000[,-1]
MATRIX_DAY0_9000_fit = MATRIX_DAY0_9000[,-1]
MATRIX_DAY0_10000_fit = MATRIX_DAY0_10000[,-1]

MATRIX_DAY0_to10000 = (MATRIX_DAY0_1000_fit +
                       MATRIX_DAY0_2000_fit +
                       MATRIX_DAY0_3000_fit +
                       MATRIX_DAY0_4000_fit +
                       MATRIX_DAY0_5000_fit +
                       MATRIX_DAY0_6000_fit +
                       MATRIX_DAY0_7000_fit +
                       MATRIX_DAY0_8000_fit +
                       MATRIX_DAY0_9000_fit +
                       MATRIX_DAY0_10000_fit)MATRIX_DAY0_to10000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_to10000.txt", delim = "\t",show_col_types = FALSE))

MATRIX_DAY0_11000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_11000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY0_12000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_12000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY0_13000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_13000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY0_14000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_14000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY0_15000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_15000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY0_16000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_16000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY0_17000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_17000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY0_18000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_18000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY0_19000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_19000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY0_20000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_20000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY0_21000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_21000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY0_22000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_22000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY0_23000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_23000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY0_24000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_24000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY0_fin = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_fin_header.txt", delim = "\t",show_col_types = FALSE))

rowname = MATRIX_DAY0_11000[,1]

MATRIX_DAY0_11000_fit = MATRIX_DAY0_11000[,-1]
MATRIX_DAY0_12000_fit = MATRIX_DAY0_12000[,-1]
MATRIX_DAY0_13000_fit = MATRIX_DAY0_13000[,-1]
MATRIX_DAY0_14000_fit = MATRIX_DAY0_14000[,-1]
MATRIX_DAY0_15000_fit = MATRIX_DAY0_15000[,-1]
MATRIX_DAY0_16000_fit = MATRIX_DAY0_16000[,-1]
MATRIX_DAY0_17000_fit = MATRIX_DAY0_17000[,-1]
MATRIX_DAY0_18000_fit = MATRIX_DAY0_18000[,-1]
MATRIX_DAY0_19000_fit = MATRIX_DAY0_19000[,-1]
MATRIX_DAY0_20000_fit = MATRIX_DAY0_20000[,-1]
MATRIX_DAY0_21000_fit = MATRIX_DAY0_21000[,-1]
MATRIX_DAY0_22000_fit = MATRIX_DAY0_22000[,-1]
MATRIX_DAY0_23000_fit = MATRIX_DAY0_23000[,-1]
MATRIX_DAY0_24000_fit = MATRIX_DAY0_24000[,-1]
MATRIX_DAY0_fin_fit = MATRIX_DAY0_fin[,-1]

MATRIX_DAY0 = (MATRIX_DAY0_to10000 +
                       MATRIX_DAY0_11000_fit +
                       MATRIX_DAY0_12000_fit +
                       MATRIX_DAY0_13000_fit +
                       MATRIX_DAY0_14000_fit +
                       MATRIX_DAY0_15000_fit +
                       MATRIX_DAY0_16000_fit +
                       MATRIX_DAY0_17000_fit +
                       MATRIX_DAY0_18000_fit +
                       MATRIX_DAY0_19000_fit +
                       MATRIX_DAY0_20000_fit +
                       MATRIX_DAY0_21000_fit +
                       MATRIX_DAY0_22000_fit +
                       MATRIX_DAY0_23000_fit +
                       MATRIX_DAY0_24000_fit +
                       MATRIX_DAY0_fin_fit)

write.table(MATRIX_DAY0,
            file = "/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_UP.txt",
            col.names=T, row.names=F, quote=F, sep="\t")

MATRIX_DAY0_UP = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_UP.txt", delim = "\t",show_col_types = FALSE))

MATRIX_DAY0_UP_MAT <- as.matrix(MATRIX_DAY0_UP)  
MATRIX_DAY0 <- MATRIX_DAY0_UP_MAT
MATRIX_DAY0[lower.tri(MATRIX_DAY0)] = t(MATRIX_DAY0)[lower.tri(MATRIX_DAY0)]
MATRIX_DAY0 = as.data.frame(MATRIX_DAY0)
write.table(MATRIX_DAY0,
            file = "/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0.txt",
            col.names=T, row.names=F, quote=F, sep="\t")MATRIX_DAY0 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0.txt", delim = "\t",show_col_types = FALSE))
rownames(MATRIX_DAY0) = c(1:nrow(MATRIX_DAY0))

MATRIX_DAY0_1000_header = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_1000_header.txt", delim = "\t",show_col_types = FALSE))
header = MATRIX_DAY0_1000_header[,1]
header = as.data.frame(header)
colnames(header) = "bin"

MATRIX_DAY0 = cbind(MATRIX_DAY0,header)

BIN_AB_data_NOREts1_p300 = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DATA/BIN_AB_data_NOREts1_p300.txt",sep = "\t", header=F)

BIN_AB_data_NOREts1_p300$bin = paste(BIN_AB_data_NOREts1_p300[,1],BIN_AB_data_NOREts1_p300[,2],BIN_AB_data_NOREts1_p300[,3])
BIN_AB_data_NOREts1_p300$Ets1 = 1

BIN_AB_data_NOREts1_p300 = BIN_AB_data_NOREts1_p300[,c(12,13,4:11)]


MATRIX_DAY0_NOR = full_join(MATRIX_DAY0,BIN_AB_data_NOREts1_p300,by="bin")
MATRIX_DAY0_NOR[is.na(MATRIX_DAY0_NOR)] = 0

Ets1_num = rownames(MATRIX_DAY0_NOR[which(MATRIX_DAY0_NOR$Ets1 == "1"),])

Ets1_num_int = as.data.frame(Ets1_num)
Ets1_num_int[,1] = as.integer(Ets1_num_int[,1])

MATRIX_DAY0_NOR_filltered = MATRIX_DAY0_NOR[Ets1_num_int[,1],]

MATRIX_DAY0_NOR_filltered1 = MATRIX_DAY0_NOR_filltered[,Ets1_num_int[,1]]
MATRIX_DAY0_NOR_filltered2 = MATRIX_DAY0_NOR_filltered[,c(24640:24649)]

MATRIX_DAY0_NOR_filltered = cbind(MATRIX_DAY0_NOR_filltered1,
                                  MATRIX_DAY0_NOR_filltered2)

MAT_forUMAP_DAY0 = MATRIX_DAY0_NOR_filltered[1:nrow(MATRIX_DAY0_NOR_filltered),1:nrow(MATRIX_DAY0_NOR_filltered)]

write.table(MAT_forUMAP_DAY0,
            file = "/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MAT_forUMAP_DAY0.txt",
            col.names=T, row.names=F, quote=F, sep="\t")
write.table(MATRIX_DAY0_NOR_filltered2,
            file = "/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0_NOR_filltered2.txt",
            col.names=T, row.names=F, quote=F, sep="\t")Rscript UMAP_function.R "euclidean" 15 2 0.01 1 15
Rscript UMAP_function.R "cosine" 15 2 0.01 1 15
Rscript UMAP_function.R "manhattan" 15 2 0.01 1 15
Rscript UMAP_function.R "hamming" 15 2 0.01 1 15
Rscript UMAP_function.R "correlation" 15 2 0.01 1 15
Rscript UMAP_function.R "categorical" 15 2 0.01 1 15MATRIX_DAY4_2000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY4_2000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY4_4000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY4_4000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY4_6000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY4_6000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY4_8000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY4_8000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY4_10000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY4_10000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY4_12000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY4_12000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY4_14000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY4_14000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY4_16000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY4_16000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY4_18000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY4_18000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY4_20000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY4_20000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY4_22000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY4_22000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY4_24000 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY4_24000_header.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY4_fin = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY4_fin_header.txt", delim = "\t",show_col_types = FALSE))

rowname = MATRIX_DAY4_2000[,1]

MATRIX_DAY4_2000_fit = MATRIX_DAY4_2000[,-1]
MATRIX_DAY4_4000_fit = MATRIX_DAY4_4000[,-1]
MATRIX_DAY4_6000_fit = MATRIX_DAY4_6000[,-1]
MATRIX_DAY4_8000_fit = MATRIX_DAY4_8000[,-1]
MATRIX_DAY4_10000_fit = MATRIX_DAY4_10000[,-1]
MATRIX_DAY4_12000_fit = MATRIX_DAY4_12000[,-1]
MATRIX_DAY4_14000_fit = MATRIX_DAY4_14000[,-1]
MATRIX_DAY4_16000_fit = MATRIX_DAY4_16000[,-1]
MATRIX_DAY4_18000_fit = MATRIX_DAY4_18000[,-1]
MATRIX_DAY4_20000_fit = MATRIX_DAY4_20000[,-1]
MATRIX_DAY4_22000_fit = MATRIX_DAY4_22000[,-1]
MATRIX_DAY4_24000_fit = MATRIX_DAY4_24000[,-1]
MATRIX_DAY4_fin_fit = MATRIX_DAY4_fin[,-1]

MATRIX_DAY4 = (MATRIX_DAY4_2000_fit +
               MATRIX_DAY4_4000_fit +
               MATRIX_DAY4_6000_fit +
               MATRIX_DAY4_8000_fit +
               MATRIX_DAY4_10000_fit +
               MATRIX_DAY4_12000_fit +
               MATRIX_DAY4_14000_fit +
               MATRIX_DAY4_16000_fit +
               MATRIX_DAY4_18000_fit +
               MATRIX_DAY4_20000_fit +
               MATRIX_DAY4_22000_fit +
               MATRIX_DAY4_24000_fit +
                       MATRIX_DAY4_fin_fit)

write.table(MATRIX_DAY4,
            file = "/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY4_UP.txt",
            col.names=T, row.names=F, quote=F, sep="\t")

MATRIX_DAY4_UP = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY4_UP.txt", delim = "\t",show_col_types = FALSE))

MATRIX_DAY4_UP_MAT <- as.matrix(MATRIX_DAY4_UP)  
MATRIX_DAY4 <- MATRIX_DAY4_UP_MAT
MATRIX_DAY4[lower.tri(MATRIX_DAY4)] = t(MATRIX_DAY4)[lower.tri(MATRIX_DAY4)]
MATRIX_DAY4 = as.data.frame(MATRIX_DAY4)
write.table(MATRIX_DAY4,
            file = "/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY4.txt",
            col.names=T, row.names=F, quote=F, sep="\t")MATRIX_DAY0 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0.txt", delim = "\t",show_col_types = FALSE))
MATRIX_DAY4 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY4.txt", delim = "\t",show_col_types = FALSE))

MATRIX_DIFF = MATRIX_DAY4 - MATRIX_DAY0

write.table(MATRIX_DIFF,
            file = "/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DIFF.txt",
            col.names=T, row.names=F, quote=F, sep="\t")# bedtools intersect -a BIN_AB_data_tmp.txt -b DAY0_Ets1_151_fit.bedpe -wa -wb > BIN_AB_Ets1_data_tmp.txt

BIN_AB_data_tmp = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DATA/full/BIN_AB_data_tmp.txt",sep = "\t", header=F)
BIN_AB_Ets1_data_tmp = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DATA/full/BIN_AB_Ets1_data_tmp.txt",sep = "\t", header=F)

BIN_AB_Ets1_data = BIN_AB_data_tmp[,1:7]
BIN_AB_Ets1_data = unique(BIN_AB_Ets1_data)

BIN_AB_Ets1_data$label = paste(BIN_AB_Ets1_data[,1],BIN_AB_Ets1_data[,2],BIN_AB_Ets1_data[,3])
BIN_AB_Ets1_data_tmp$label = paste(BIN_AB_Ets1_data_tmp[,1],BIN_AB_Ets1_data_tmp[,2],BIN_AB_Ets1_data_tmp[,3])

BIN_AB_Ets1_data$score = 0

for( i in 1:nrow(BIN_AB_Ets1_data)){
BIN_AB_Ets1_data[i,9] = nrow(BIN_AB_Ets1_data_tmp[which(BIN_AB_Ets1_data_tmp$label == BIN_AB_Ets1_data[i,8]),])
    }

BIN_AB_Ets1_data = BIN_AB_Ets1_data[order(BIN_AB_Ets1_data$score),]

BIN_AB_Ets1_data$Q = "NOT"
GAP = 2464

for(i in 1:9){
BIN_AB_Ets1_data[((GAP*(i-1))+1):(GAP*i),10] = paste("Q",i,sep="")
    }
BIN_AB_Ets1_data[((2464*9)+1):nrow(BIN_AB_Ets1_data),10] = "Q10"

write.table(BIN_AB_Ets1_data,
            file = "/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DATA/full/BIN_AB_Ets1_data.txt",
            col.names=F, row.names=F, quote=F, sep="\t")# bedtools intersect -a BIN_AB_Ets1_data.txt -b DAY0_p300_151_fit.bedpe -wa -wb > BIN_AB_Ets1_p300_data_tmp.txt

BIN_AB_Ets1_data = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DATA/full/BIN_AB_Ets1_data.txt",sep = "\t", header=F)
BIN_AB_Ets1_p300_data_tmp = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DATA/full/BIN_AB_Ets1_p300_data_tmp.txt",sep = "\t", header=F)

BIN_AB_Ets1_data = unique(BIN_AB_Ets1_data)
BIN_AB_Ets1_data$score = 0

for( i in 1:nrow(BIN_AB_Ets1_data)){
BIN_AB_Ets1_data[i,11] = nrow(BIN_AB_Ets1_p300_data_tmp[which(BIN_AB_Ets1_p300_data_tmp[,8] == BIN_AB_Ets1_data[i,8]),])
    }

BIN_AB_Ets1_data = BIN_AB_Ets1_data[order(BIN_AB_Ets1_data$score),]

BIN_AB_Ets1_data$Q = "NOT"
GAP = 2464

for(i in 1:9){
BIN_AB_Ets1_data[((GAP*(i-1))+1):(GAP*i),12] = paste("Q",i,sep="")
    }
BIN_AB_Ets1_data[((2464*9)+1):nrow(BIN_AB_Ets1_data),12] = "Q10"

BIN_AB_Ets1_data = BIN_AB_Ets1_data[,-8]

write.table(BIN_AB_Ets1_data,
            file = "/data3/psg/NGS_2025/4DN/HiC_data/FiG5/DATA/full/BIN_AB_Ets1_p300_data.txt",
            col.names=F, row.names=F, quote=F, sep="\t")

```R

```


```R

```


```R
options(repr.plot.width = 5, repr.plot.height = 5, repr.plot.res = 1000, repr.plot.pointsize = 40)

MAP = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY0.txt", delim = "\t",show_col_types = FALSE))

conte=MAP
conte = conte +1
conte = log2(conte)
rownames(conte)=c(1:nrow(MAP))

PCA=conte
fig = pheatmap(PCA,
               cluster_rows = F,
               cluster_cols = F,
               show_rownames = F,
               show_colnames = F,
               legend = FALSE,
               border_color = FALSE,
               color = colorRampPalette(c("white","red","red","darkred"))(100),
               breaks=seq(from=0,to=10,length.out = 100))
```


    
![png](output_17_0.png)
    



```R
options(repr.plot.width = 5, repr.plot.height = 5, repr.plot.res = 1000, repr.plot.pointsize = 40)

MAP = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DAY4.txt", delim = "\t",show_col_types = FALSE))

conte=MAP
conte = conte +1
conte = log2(conte)
rownames(conte)=c(1:nrow(MAP))

PCA=conte
fig = pheatmap(PCA,
               cluster_rows = F,
               cluster_cols = F,
               show_rownames = F,
               show_colnames = F,
               legend = FALSE,
               border_color = FALSE,
               color = colorRampPalette(c("white","red","red","darkred"))(100),
               breaks=seq(from=0,to=10,length.out = 100))
```


    
![png](output_18_0.png)
    



```R
options(repr.plot.width = 5, repr.plot.height = 5, repr.plot.res = 1000, repr.plot.pointsize = 40)

MAP = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/MATRIX/MATRIX_DIFF.txt", delim = "\t",show_col_types = FALSE))

conte=MAP
#conte = conte +1
#conte = log2(conte)
rownames(conte)=c(1:nrow(MAP))

PCA=conte
fig = pheatmap(PCA,
               cluster_rows = F,
               cluster_cols = F,
               show_rownames = F,
               show_colnames = F,
               legend = FALSE,
               border_color = FALSE,
               color = colorRampPalette(c("darkblue","white","darkred"))(100),
               breaks=seq(from=-3,to=3,length.out = 100))
```


    
![png](output_19_0.png)
    



```R
ppp=arrangeGrob(

UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_50_2_01_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_cosine_50_2_01_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_manhattan_50_2_01_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_hamming_50_2_01_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_50_2_01_1.txt"),

    ncol = 5,
    nrow = 1)

options(repr.plot.width = 40, repr.plot.height = 8, repr.plot.res = 200, repr.plot.pointsize = 40)
grid.arrange(ppp)
```


    
![png](output_20_0.png)
    



```R
ppp=arrangeGrob(

UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_500_2_00001_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_100_2_00001_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_50_2_00001_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_10_2_00001_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_5_2_00001_1.txt"),
    
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_500_2_0001_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_100_2_0001_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_50_2_0001_1.txt"), 
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_10_2_0001_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_5_2_0001_1.txt"),

UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_500_2_001_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_100_2_001_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_50_2_001_1.txt"), 
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_10_2_001_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_5_2_001_1.txt"),

UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_500_2_01_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_100_2_01_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_50_2_01_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_10_2_01_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_5_2_01_1.txt"),

UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_500_2_1_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_100_2_1_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_50_2_1_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_10_2_1_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_5_2_1_1.txt"),

    ncol = 5,
    nrow = 5)

options(repr.plot.width = 40, repr.plot.height = 40, repr.plot.res = 200, repr.plot.pointsize = 40)
grid.arrange(ppp)
```


    
![png](output_21_0.png)
    



```R
ppp=arrangeGrob(

UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_500_2_00001_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_100_2_00001_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_50_2_00001_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_10_2_00001_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_5_2_00001_1.txt"),
    
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_500_2_0001_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_100_2_0001_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_50_2_0001_1.txt"), 
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_10_2_0001_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_5_2_0001_1.txt"),

UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_500_2_001_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_100_2_001_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_50_2_001_1.txt"), 
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_10_2_001_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_5_2_001_1.txt"),

UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_500_2_01_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_100_2_01_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_50_2_01_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_10_2_01_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_5_2_01_1.txt"),

UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_500_2_1_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_100_2_1_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_50_2_1_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_10_2_1_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_5_2_1_1.txt"),

    ncol = 5,
    nrow = 5)

options(repr.plot.width = 40, repr.plot.height = 40, repr.plot.res = 200, repr.plot.pointsize = 40)
grid.arrange(ppp)
```


    
![png](output_22_0.png)
    



```R
ppp=arrangeGrob(

UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_10_2_1_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_10_2_1_3.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_10_2_1_5.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_10_2_1_10.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_10_2_1_50.txt"),


    ncol = 5,
    nrow = 1)

options(repr.plot.width = 40, repr.plot.height = 8, repr.plot.res = 200, repr.plot.pointsize = 40)
grid.arrange(ppp)
```


    
![png](output_23_0.png)
    



```R
ppp=arrangeGrob(

UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_50_2_001_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_50_2_001_3.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_50_2_001_5.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_50_2_001_10.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_50_2_001_50.txt"),


    ncol = 5,
    nrow = 1)

options(repr.plot.width = 40, repr.plot.height = 8, repr.plot.res = 200, repr.plot.pointsize = 40)
grid.arrange(ppp)
```


    
![png](output_24_0.png)
    



```R
ppp=arrangeGrob(

UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_50_2_1_1.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_50_2_1_3.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_50_2_1_5.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_50_2_1_10.txt"),
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_50_2_1_50.txt"),


    ncol = 5,
    nrow = 1)

options(repr.plot.width = 40, repr.plot.height = 8, repr.plot.res = 200, repr.plot.pointsize = 40)
grid.arrange(ppp)
```


    
![png](output_25_0.png)
    



```R
options(repr.plot.width = 10, repr.plot.height = 10, repr.plot.res = 1000, repr.plot.pointsize = 40)
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_50_2_1_1.txt")
UMAP_viewer("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_euclidean_10_2_1_1.txt")
```


    
![png](output_26_0.png)
    



    
![png](output_26_1.png)
    



```R
set.seed(7777)

UMAP_DIFF_correlation_50_2_1_1 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_50_2_1_1.txt", delim = "\t",show_col_types = FALSE))
UMAP_DIFF_correlation_50_2_1_1_fit = UMAP_DIFF_correlation_50_2_1_1[,1:2]
UMAP_DIFF_correlation_50_2_1_1_label = paste(UMAP_DIFF_correlation_50_2_1_1[,3],UMAP_DIFF_correlation_50_2_1_1[,4],UMAP_DIFF_correlation_50_2_1_1[,5])
UMAP_DIFF_correlation_50_2_1_1_label = as.data.frame(UMAP_DIFF_correlation_50_2_1_1_label)
colnames(UMAP_DIFF_correlation_50_2_1_1_label) = "label"

result <- kmeans(UMAP_DIFF_correlation_50_2_1_1_fit, centers = 22, iter.max = 10000, algorithm = "MacQueen")

UMAP_DIFF_correlation_50_2_1_1 = cbind(UMAP_DIFF_correlation_50_2_1_1,as.data.frame(result$cluster))
colnames(UMAP_DIFF_correlation_50_2_1_1)[ncol(UMAP_DIFF_correlation_50_2_1_1)] = "cluster"

i = 1

ABsignal_cluster = cbind(i,median(UMAP_DIFF_correlation_50_2_1_1[which(UMAP_DIFF_correlation_50_2_1_1$cluster == i),8]))
ABsignal_cluster = as.data.frame(ABsignal_cluster)
colnames(ABsignal_cluster) = c("cluster","ABsignal")

for(i in 2:22){
    ABsignal_cluster_frag = cbind(i,median(UMAP_DIFF_correlation_50_2_1_1[which(UMAP_DIFF_correlation_50_2_1_1$cluster == i),8]))
    ABsignal_cluster_frag = as.data.frame(ABsignal_cluster_frag)
    colnames(ABsignal_cluster_frag) = c("cluster","ABsignal")

    ABsignal_cluster = rbind(ABsignal_cluster,ABsignal_cluster_frag)
    }

ABsignal_cluster = ABsignal_cluster[order(-ABsignal_cluster$ABsignal),]
rownames(ABsignal_cluster) = NULL
ABsignal_cluster[,1] = as.factor(ABsignal_cluster[,1])

UMAP_DIFF_correlation_50_2_1_1_forPlot_DIFF = UMAP_DIFF_correlation_50_2_1_1[,c(14,8)]

UMAP_DIFF_correlation_50_2_1_1_forPlot_Ets1 = UMAP_DIFF_correlation_50_2_1_1[,c(14,10)]

UMAP_DIFF_correlation_50_2_1_1_forPlot_Ets1[,2] = log2(UMAP_DIFF_correlation_50_2_1_1_forPlot_Ets1[,2]+1)

head(UMAP_DIFF_correlation_50_2_1_1)
```


<table class="dataframe">
<caption>A data.frame: 6 × 14</caption>
<thead>
	<tr><th></th><th scope=col>X</th><th scope=col>Y</th><th scope=col>chr</th><th scope=col>str</th><th scope=col>end</th><th scope=col>AB_signal_DAY0</th><th scope=col>AB_signal_DAY4</th><th scope=col>AB_signal_DIFF</th><th scope=col>DIFF_Q</th><th scope=col>Ets1</th><th scope=col>Ets1_Q</th><th scope=col>p300</th><th scope=col>p300_Q</th><th scope=col>cluster</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>17.71553</td><td>-3.816955</td><td>chr1</td><td>     1</td><td>100000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><th scope=row>2</th><td>17.44521</td><td>-3.702754</td><td>chr1</td><td>100001</td><td>200000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><th scope=row>3</th><td>18.01741</td><td>-3.566366</td><td>chr1</td><td>200001</td><td>300000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><th scope=row>4</th><td>17.51338</td><td>-3.936230</td><td>chr1</td><td>300001</td><td>400000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><th scope=row>5</th><td>16.92375</td><td>-3.072998</td><td>chr1</td><td>400001</td><td>500000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><th scope=row>6</th><td>17.48633</td><td>-3.129353</td><td>chr1</td><td>500001</td><td>600000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
</tbody>
</table>




```R
UMAP_DIFF_correlation_50_2_1_1
```


<table class="dataframe">
<caption>A data.frame: 24639 × 14</caption>
<thead>
	<tr><th scope=col>X</th><th scope=col>Y</th><th scope=col>chr</th><th scope=col>str</th><th scope=col>end</th><th scope=col>AB_signal_DAY0</th><th scope=col>AB_signal_DAY4</th><th scope=col>AB_signal_DIFF</th><th scope=col>DIFF_Q</th><th scope=col>Ets1</th><th scope=col>Ets1_Q</th><th scope=col>p300</th><th scope=col>p300_Q</th><th scope=col>cluster</th></tr>
	<tr><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>17.71553</td><td>-3.816955</td><td>chr1</td><td>      1</td><td> 100000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>17.44521</td><td>-3.702754</td><td>chr1</td><td> 100001</td><td> 200000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>18.01741</td><td>-3.566366</td><td>chr1</td><td> 200001</td><td> 300000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>17.51338</td><td>-3.936230</td><td>chr1</td><td> 300001</td><td> 400000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>16.92375</td><td>-3.072998</td><td>chr1</td><td> 400001</td><td> 500000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>17.48633</td><td>-3.129353</td><td>chr1</td><td> 500001</td><td> 600000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>17.91799</td><td>-3.588283</td><td>chr1</td><td> 600001</td><td> 700000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>17.63700</td><td>-4.019646</td><td>chr1</td><td> 700001</td><td> 800000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>17.38618</td><td>-3.690047</td><td>chr1</td><td> 800001</td><td> 900000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>18.27905</td><td>-3.380016</td><td>chr1</td><td> 900001</td><td>1000000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>17.55659</td><td>-3.587095</td><td>chr1</td><td>1000001</td><td>1100000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>17.39846</td><td>-4.293754</td><td>chr1</td><td>1100001</td><td>1200000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>17.88053</td><td>-3.458113</td><td>chr1</td><td>1200001</td><td>1300000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>17.82925</td><td>-3.909251</td><td>chr1</td><td>1300001</td><td>1400000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>17.92680</td><td>-3.034559</td><td>chr1</td><td>1400001</td><td>1500000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>18.47013</td><td>-2.995713</td><td>chr1</td><td>1500001</td><td>1600000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>18.11814</td><td>-2.738005</td><td>chr1</td><td>1600001</td><td>1700000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>17.70843</td><td>-3.734890</td><td>chr1</td><td>1700001</td><td>1800000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>18.17000</td><td>-3.351405</td><td>chr1</td><td>1800001</td><td>1900000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>17.31908</td><td>-4.013596</td><td>chr1</td><td>1900001</td><td>2000000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>17.74782</td><td>-2.661716</td><td>chr1</td><td>2000001</td><td>2100000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>17.41088</td><td>-2.745890</td><td>chr1</td><td>2100001</td><td>2200000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>18.38206</td><td>-3.047066</td><td>chr1</td><td>2200001</td><td>2300000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>18.37890</td><td>-3.419734</td><td>chr1</td><td>2300001</td><td>2400000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>17.82785</td><td>-3.701291</td><td>chr1</td><td>2400001</td><td>2500000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>17.95631</td><td>-3.197901</td><td>chr1</td><td>2500001</td><td>2600000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>18.26710</td><td>-2.786247</td><td>chr1</td><td>2600001</td><td>2700000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>17.52663</td><td>-3.777373</td><td>chr1</td><td>2700001</td><td>2800000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>17.69055</td><td>-3.871503</td><td>chr1</td><td>2800001</td><td>2900000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>16.69961</td><td>-3.629432</td><td>chr1</td><td>2900001</td><td>3000000</td><td>0.0005551651</td><td>0.0002517516</td><td>-0.0003034136</td><td>Q6</td><td>0</td><td>Q1</td><td>0</td><td>Q1</td><td>1</td></tr>
	<tr><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><td> 3.186826</td><td> 2.555259</td><td>chr9</td><td>121600001</td><td>121700000</td><td> 0.012813151</td><td> 0.016432221807</td><td> 0.0036190706</td><td>Q7</td><td>59</td><td>Q1</td><td>58</td><td>Q3</td><td>12</td></tr>
	<tr><td> 3.575922</td><td> 2.922670</td><td>chr9</td><td>121700001</td><td>121800000</td><td> 0.012813151</td><td> 0.022693115272</td><td> 0.0098799641</td><td>Q8</td><td>59</td><td>Q1</td><td>48</td><td>Q1</td><td>12</td></tr>
	<tr><td> 3.724544</td><td> 2.911502</td><td>chr9</td><td>121800001</td><td>121900000</td><td> 0.012813151</td><td> 0.018328654528</td><td> 0.0055155033</td><td>Q7</td><td>45</td><td>Q1</td><td>60</td><td>Q3</td><td>12</td></tr>
	<tr><td> 3.444619</td><td> 2.710122</td><td>chr9</td><td>121900001</td><td>122000000</td><td> 0.012813151</td><td> 0.015798759741</td><td> 0.0029856085</td><td>Q7</td><td>60</td><td>Q1</td><td>47</td><td>Q1</td><td>12</td></tr>
	<tr><td> 3.292691</td><td> 2.593361</td><td>chr9</td><td>122000001</td><td>122100000</td><td> 0.012813151</td><td> 0.011423359895</td><td>-0.0013897913</td><td>Q5</td><td>60</td><td>Q1</td><td>56</td><td>Q3</td><td>12</td></tr>
	<tr><td> 3.655697</td><td> 2.522079</td><td>chr9</td><td>122100001</td><td>122200000</td><td> 0.012813151</td><td> 0.013355886355</td><td> 0.0005427352</td><td>Q6</td><td>74</td><td>Q4</td><td>54</td><td>Q2</td><td>12</td></tr>
	<tr><td> 3.461131</td><td> 2.912646</td><td>chr9</td><td>122200001</td><td>122300000</td><td> 0.012813151</td><td> 0.013718590809</td><td> 0.0009054396</td><td>Q6</td><td>81</td><td>Q6</td><td>60</td><td>Q4</td><td>12</td></tr>
	<tr><td> 3.449029</td><td> 2.658868</td><td>chr9</td><td>122300001</td><td>122400000</td><td> 0.012813151</td><td> 0.014753547649</td><td> 0.0019403964</td><td>Q7</td><td>78</td><td>Q5</td><td>72</td><td>Q7</td><td>12</td></tr>
	<tr><td> 3.677232</td><td> 2.895640</td><td>chr9</td><td>122400001</td><td>122500000</td><td> 0.012813151</td><td> 0.000645586873</td><td>-0.0121675643</td><td>Q3</td><td>72</td><td>Q3</td><td>79</td><td>Q8</td><td>12</td></tr>
	<tr><td> 3.618804</td><td> 2.955827</td><td>chr9</td><td>122500001</td><td>122600000</td><td>-0.002814566</td><td>-0.000007013667</td><td> 0.0028075526</td><td>Q7</td><td>67</td><td>Q2</td><td>49</td><td>Q1</td><td>12</td></tr>
	<tr><td> 3.451022</td><td> 3.003456</td><td>chr9</td><td>122600001</td><td>122700000</td><td> 0.008611739</td><td>-0.001481119049</td><td>-0.0100928581</td><td>Q3</td><td>77</td><td>Q5</td><td>54</td><td>Q2</td><td>12</td></tr>
	<tr><td> 3.951870</td><td> 2.812284</td><td>chr9</td><td>122700001</td><td>122800000</td><td> 0.008611739</td><td> 0.002196627116</td><td>-0.0064151119</td><td>Q4</td><td>71</td><td>Q3</td><td>56</td><td>Q3</td><td>12</td></tr>
	<tr><td> 3.432698</td><td> 3.151766</td><td>chr9</td><td>122800001</td><td>122900000</td><td> 0.008611739</td><td> 0.009296500027</td><td> 0.0006847610</td><td>Q6</td><td>61</td><td>Q2</td><td>52</td><td>Q2</td><td>12</td></tr>
	<tr><td> 3.266199</td><td> 2.931793</td><td>chr9</td><td>122900001</td><td>123000000</td><td> 0.008611739</td><td> 0.007433401331</td><td>-0.0011783377</td><td>Q5</td><td>77</td><td>Q5</td><td>48</td><td>Q1</td><td>12</td></tr>
	<tr><td> 3.342125</td><td> 3.046007</td><td>chr9</td><td>123000001</td><td>123100000</td><td> 0.008611739</td><td> 0.011846029719</td><td> 0.0032342907</td><td>Q7</td><td>67</td><td>Q2</td><td>57</td><td>Q3</td><td>12</td></tr>
	<tr><td> 3.554562</td><td> 3.151502</td><td>chr9</td><td>123100001</td><td>123200000</td><td> 0.008611739</td><td> 0.011835029860</td><td> 0.0032232909</td><td>Q7</td><td>75</td><td>Q4</td><td>55</td><td>Q2</td><td>12</td></tr>
	<tr><td> 3.710124</td><td> 3.036774</td><td>chr9</td><td>123200001</td><td>123300000</td><td> 0.008611739</td><td> 0.008025082371</td><td>-0.0005866566</td><td>Q6</td><td>64</td><td>Q2</td><td>44</td><td>Q1</td><td>12</td></tr>
	<tr><td> 3.332471</td><td> 2.862528</td><td>chr9</td><td>123300001</td><td>123400000</td><td> 0.008611739</td><td> 0.006959073629</td><td>-0.0016526654</td><td>Q5</td><td>67</td><td>Q2</td><td>52</td><td>Q2</td><td>12</td></tr>
	<tr><td> 3.263835</td><td> 3.216507</td><td>chr9</td><td>123400001</td><td>123500000</td><td> 0.008611739</td><td> 0.013373870249</td><td> 0.0047621312</td><td>Q7</td><td>58</td><td>Q1</td><td>53</td><td>Q2</td><td>12</td></tr>
	<tr><td> 3.263149</td><td> 3.013936</td><td>chr9</td><td>123500001</td><td>123600000</td><td> 0.008611739</td><td> 0.017241396668</td><td> 0.0086296577</td><td>Q8</td><td>65</td><td>Q2</td><td>55</td><td>Q2</td><td>12</td></tr>
	<tr><td> 3.184123</td><td> 2.363146</td><td>chr9</td><td>123600001</td><td>123700000</td><td> 0.008611739</td><td> 0.008357317251</td><td>-0.0002544218</td><td>Q6</td><td>80</td><td>Q5</td><td>52</td><td>Q2</td><td>12</td></tr>
	<tr><td> 3.108238</td><td> 2.581162</td><td>chr9</td><td>123700001</td><td>123800000</td><td> 0.008611739</td><td> 0.009516738774</td><td> 0.0009049998</td><td>Q6</td><td>68</td><td>Q3</td><td>76</td><td>Q8</td><td>12</td></tr>
	<tr><td> 2.878594</td><td> 2.484939</td><td>chr9</td><td>123800001</td><td>123900000</td><td> 0.008611739</td><td> 0.006162917629</td><td>-0.0024488214</td><td>Q5</td><td>49</td><td>Q1</td><td>37</td><td>Q1</td><td>12</td></tr>
	<tr><td> 2.999261</td><td> 2.783317</td><td>chr9</td><td>123900001</td><td>124000000</td><td>-0.023115649</td><td>-0.019120566773</td><td> 0.0039950823</td><td>Q7</td><td>68</td><td>Q3</td><td>40</td><td>Q1</td><td>12</td></tr>
	<tr><td> 3.163915</td><td> 2.493255</td><td>chr9</td><td>124000001</td><td>124100000</td><td>-0.023115649</td><td>-0.028456103551</td><td>-0.0053404545</td><td>Q4</td><td>67</td><td>Q2</td><td>59</td><td>Q3</td><td>12</td></tr>
	<tr><td> 3.272448</td><td> 2.143018</td><td>chr9</td><td>124100001</td><td>124200000</td><td>-0.023115649</td><td>-0.010599677910</td><td> 0.0125159712</td><td>Q8</td><td>44</td><td>Q1</td><td>35</td><td>Q1</td><td>12</td></tr>
	<tr><td> 3.016973</td><td> 2.320808</td><td>chr9</td><td>124200001</td><td>124300000</td><td>-0.023115649</td><td>-0.010114098027</td><td> 0.0130015510</td><td>Q8</td><td>17</td><td>Q1</td><td>23</td><td>Q1</td><td>12</td></tr>
	<tr><td> 2.795636</td><td> 2.273545</td><td>chr9</td><td>124300001</td><td>124400000</td><td>-0.023115649</td><td>-0.016978722687</td><td> 0.0061369264</td><td>Q7</td><td>53</td><td>Q1</td><td>40</td><td>Q1</td><td>12</td></tr>
	<tr><td> 3.443125</td><td> 2.576841</td><td>chr9</td><td>124400001</td><td>124500000</td><td>-0.023115649</td><td>-0.005525618196</td><td> 0.0175900309</td><td>Q9</td><td>55</td><td>Q1</td><td>46</td><td>Q1</td><td>12</td></tr>
	<tr><td>17.829579</td><td>-3.718128</td><td>chr9</td><td>124500001</td><td>124600000</td><td> 0.000000000</td><td> 0.000000000000</td><td> 0.0000000000</td><td>Q6</td><td> 0</td><td>Q1</td><td> 0</td><td>Q1</td><td> 1</td></tr>
</tbody>
</table>




```R
options(repr.plot.width = 15, repr.plot.height = 7, repr.plot.res = 1000, repr.plot.pointsize = 40)

ggplot(UMAP_DIFF_correlation_50_2_1_1_forPlot_DIFF, aes(x=factor(cluster), y=AB_signal_DIFF, fill=factor(cluster)))+
scale_x_discrete(limits=ABsignal_cluster[,1])+
geom_boxplot(#fill=c("dodgerblue3","yellow3"),
             size=1,
             color='black',
             width=0.5,
             alpha=0.7)+
geom_violin(adjust=1,
            scale='count',
            fill=c("black"),
            alpha=0.5)+
    guides(fill="none")+
theme_minimal()+
    labs(title="",
        subtitle = "")+
    xlab("")+
    ylab("AB compartment singal\n")+
theme(plot.title=element_text(face="bold.italic",hjust=0.5,size=25,color = "black"),
      axis.text.x=element_text(face="bold",size=12,color = "black"),
      axis.text.y=element_text(face="bold",size=10),
      plot.subtitle=element_text(face="bold",hjust=0.5,vjust=-3,size=15,color = "black"),
      axis.title.x = element_text(size=0,color = "black"),
      axis.title.y = element_text(face="bold",size=15,color = "black"))
```


    
![png](output_29_0.png)
    



```R
options(repr.plot.width = 15, repr.plot.height = 7, repr.plot.res = 1000, repr.plot.pointsize = 40)

ggplot(UMAP_DIFF_correlation_50_2_1_1_forPlot_Ets1, aes(x=factor(cluster), y=Ets1, fill=factor(cluster)))+
scale_x_discrete(limits=ABsignal_cluster[,1])+
geom_boxplot(#fill=c("dodgerblue3","yellow3"),
             size=1,
             color='black',
             width=0.5,
             alpha=0.7)+
geom_violin(adjust=1,
            scale='count',
            fill=c("black"),
            alpha=0.5)+
    guides(fill="none")+
theme_minimal()+
    labs(title="",
        subtitle = "")+
    xlab("")+
    ylab("AB compartment singal\n")+
theme(plot.title=element_text(face="bold.italic",hjust=0.5,size=25,color = "black"),
      axis.text.x=element_text(face="bold",size=12,color = "black"),
      axis.text.y=element_text(face="bold",size=10),
      plot.subtitle=element_text(face="bold",hjust=0.5,vjust=-3,size=15,color = "black"),
      axis.title.x = element_text(size=0,color = "black"),
      axis.title.y = element_text(face="bold",size=15,color = "black"))
```


    
![png](output_30_0.png)
    



```R

```


```R

```


```R

```


```R

```


```R
set.seed(7777)

UMAP_DIFF_correlation_50_2_1_1 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_50_2_1_1.txt", delim = "\t",show_col_types = FALSE))
UMAP_DIFF_correlation_50_2_1_1_fit = UMAP_DIFF_correlation_50_2_1_1[,1:2]
UMAP_DIFF_correlation_50_2_1_1_label = paste(UMAP_DIFF_correlation_50_2_1_1[,3],UMAP_DIFF_correlation_50_2_1_1[,4],UMAP_DIFF_correlation_50_2_1_1[,5])
UMAP_DIFF_correlation_50_2_1_1_label = as.data.frame(UMAP_DIFF_correlation_50_2_1_1_label)
colnames(UMAP_DIFF_correlation_50_2_1_1_label) = "label"

result <- kmeans(UMAP_DIFF_correlation_50_2_1_1_fit, centers = 22, iter.max = 10000, algorithm = "MacQueen")

UMAP_DIFF_correlation_50_2_1_1 = cbind(UMAP_DIFF_correlation_50_2_1_1,as.data.frame(result$cluster))
colnames(UMAP_DIFF_correlation_50_2_1_1)[ncol(UMAP_DIFF_correlation_50_2_1_1)] = "cluster"

options(repr.plot.width = 8, repr.plot.height = 8, repr.plot.res = 1000, repr.plot.pointsize = 30)

color_palette <- c("darkblue", "#B71C1C","#283593", "#3949AB", "#5C6BC0", "#9FA8DA",
                   "#757575", "#FFCDD2", "#E57373", "#E53935")

ggplot(UMAP_DIFF_correlation_50_2_1_1, aes(X, Y),colour = factor(DIFF_Q)) +
geom_point(data=subset(UMAP_DIFF_correlation_50_2_1_1,cluster==14), colour = "#32B2F4", size=1) + 
geom_point(data=subset(UMAP_DIFF_correlation_50_2_1_1,cluster==13), colour = "#33BADD", size=1) + 

geom_point(data=subset(UMAP_DIFF_correlation_50_2_1_1,cluster==1), colour = "grey70", size=1) +
geom_point(data=subset(UMAP_DIFF_correlation_50_2_1_1,cluster==2), colour = "grey70", size=1) +
geom_point(data=subset(UMAP_DIFF_correlation_50_2_1_1,cluster==3), colour = "grey70", size=1) +
geom_point(data=subset(UMAP_DIFF_correlation_50_2_1_1,cluster==4), colour = "grey70", size=1) +
geom_point(data=subset(UMAP_DIFF_correlation_50_2_1_1,cluster==5), colour = "grey70", size=1) +
geom_point(data=subset(UMAP_DIFF_correlation_50_2_1_1,cluster==6), colour = "grey70", size=1) +
geom_point(data=subset(UMAP_DIFF_correlation_50_2_1_1,cluster==7), colour = "grey70", size=1) +
geom_point(data=subset(UMAP_DIFF_correlation_50_2_1_1,cluster==8), colour = "grey70", size=1) +
geom_point(data=subset(UMAP_DIFF_correlation_50_2_1_1,cluster==9), colour = "grey70", size=1) +
geom_point(data=subset(UMAP_DIFF_correlation_50_2_1_1,cluster==10), colour = "grey70", size=1) +
geom_point(data=subset(UMAP_DIFF_correlation_50_2_1_1,cluster==11), colour = "grey70", size=1) +
geom_point(data=subset(UMAP_DIFF_correlation_50_2_1_1,cluster==12), colour = "grey70", size=1) +
geom_point(data=subset(UMAP_DIFF_correlation_50_2_1_1,cluster==15), colour = "grey70", size=1) +
geom_point(data=subset(UMAP_DIFF_correlation_50_2_1_1,cluster==16), colour = "grey70", size=1) +
geom_point(data=subset(UMAP_DIFF_correlation_50_2_1_1,cluster==17), colour = "grey70", size=1) +
geom_point(data=subset(UMAP_DIFF_correlation_50_2_1_1,cluster==18), colour = "grey70", size=1) +
geom_point(data=subset(UMAP_DIFF_correlation_50_2_1_1,cluster==19), colour = "grey70", size=1) +
geom_point(data=subset(UMAP_DIFF_correlation_50_2_1_1,cluster==20), colour = "grey70", size=1) +
geom_point(data=subset(UMAP_DIFF_correlation_50_2_1_1,cluster==21), colour = "grey70", size=1) +
geom_point(data=subset(UMAP_DIFF_correlation_50_2_1_1,cluster==22), colour = "grey70", size=1) +
scale_color_manual(values = color_palette) +
theme_classic(base_size = 12) +
#xlim(-20,20)+
#ylim(-20,20)+
ggtitle(NULL)+
xlab("UMAP1") +
ylab("UMAP2") +
theme(plot.title=element_text(face="bold",hjust=0.5,size=25,color = "darkblue"),
      axis.text.x=element_text(face="bold",size=20),
      axis.text.y=element_text(face="bold",size=20), 
      axis.title.x = element_text(face="bold",size = 20,color = "darkblue"),
      axis.title.y = element_text(face="bold",size = 20,color = "darkblue"),
      legend.position="none", 
      legend.title=element_text(face="bold",size=15), 
      legend.text=element_text(face="bold",size=15))
```


    
![png](output_35_0.png)
    



```R
num_cluster = 1

subdata = as.data.frame(UMAP_DIFF_correlation_50_2_1_1[which(UMAP_DIFF_correlation_50_2_1_1$cluster == num_cluster),])

AB_DIFF = c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10")

AB_DIFF = as.data.frame(AB_DIFF)
colnames(AB_DIFF) = "AB_DIFF"
AB_DIFF$score = 0

for(i in 1:nrow(AB_DIFF)){
AB_DIFF[i,2] = nrow(subdata[which(subdata$DIFF_Q == AB_DIFF[i,1]),])
    }
num_cluster_value = round(sum(AB_DIFF[7:10,2])/sum(AB_DIFF[1:5,2]),2)

DIFF_clusterDATA = as.data.frame(cbind(num_cluster,num_cluster_value))


for(s in 2:22){
num_cluster = s

subdata = as.data.frame(UMAP_DIFF_correlation_50_2_1_1[which(UMAP_DIFF_correlation_50_2_1_1$cluster == num_cluster),])

AB_DIFF = c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10")

AB_DIFF = as.data.frame(AB_DIFF)
colnames(AB_DIFF) = "AB_DIFF"
AB_DIFF$score = 0

for(i in 1:nrow(AB_DIFF)){
AB_DIFF[i,2] = nrow(subdata[which(subdata$DIFF_Q == AB_DIFF[i,1]),])
    }
num_cluster_value = round(sum(AB_DIFF[7:10,2])/sum(AB_DIFF[1:5,2]),2)

DIFF_clusterDATA_frag = as.data.frame(cbind(num_cluster,num_cluster_value))

    DIFF_clusterDATA = rbind(DIFF_clusterDATA,DIFF_clusterDATA_frag)

    }

DIFF_clusterDATA = DIFF_clusterDATA[order(-DIFF_clusterDATA$num_cluster_value),]
DIFF_clustersorted = DIFF_clusterDATA[,1]
DIFF_clustersorted = as.data.frame(DIFF_clustersorted)
colnames(DIFF_clustersorted) = "cluster"
DIFF_clustersorted[,1] = as.factor(DIFF_clustersorted[,1])

DIFF_clusterDATA$sort = c(1:22)
DIFF_clusterDATA[,3] = as.factor(DIFF_clusterDATA[,3])

DIFF_clusterDATA$label = "DIFF"

head(DIFF_clusterDATA)
```


<table class="dataframe">
<caption>A data.frame: 6 × 4</caption>
<thead>
	<tr><th></th><th scope=col>num_cluster</th><th scope=col>num_cluster_value</th><th scope=col>sort</th><th scope=col>label</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>14</th><td>14</td><td>1.76</td><td>1</td><td>DIFF</td></tr>
	<tr><th scope=row>3</th><td> 3</td><td>1.52</td><td>2</td><td>DIFF</td></tr>
	<tr><th scope=row>7</th><td> 7</td><td>1.51</td><td>3</td><td>DIFF</td></tr>
	<tr><th scope=row>16</th><td>16</td><td>1.39</td><td>4</td><td>DIFF</td></tr>
	<tr><th scope=row>17</th><td>17</td><td>1.39</td><td>5</td><td>DIFF</td></tr>
	<tr><th scope=row>8</th><td> 8</td><td>0.93</td><td>6</td><td>DIFF</td></tr>
</tbody>
</table>




```R
DIFF_clusterDATA[,1]
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>14</li><li>3</li><li>7</li><li>16</li><li>17</li><li>8</li><li>15</li><li>19</li><li>12</li><li>11</li><li>21</li><li>20</li><li>6</li><li>18</li><li>5</li><li>22</li><li>1</li><li>9</li><li>13</li><li>10</li><li>4</li><li>2</li></ol>




```R
options(repr.plot.width = 5, repr.plot.height = 3, repr.plot.res = 1000, repr.plot.pointsize = 40)
ggplot(DIFF_clusterDATA,aes(sort, y = num_cluster_value, group = "DIFF")) +
geom_line(color="black", size=0.5, alpha = 1) +
geom_point(color = "black", size = 2) + 
theme_classic(base_size = 10) +
ggtitle(NULL)+
xlab(NULL) +
ylab(NULL) +
theme(plot.title=element_text(face="bold",hjust=0.5,size=25,color = "darkblue"),
      axis.text.x=element_text(face="bold",size=0),
      axis.text.y=element_text(face="bold",size=10),
      axis.title.x = element_text(face="bold",size = 20),
      axis.title.y = element_text(face="bold",size = 15,color = "darkblue"),
      legend.title=element_text(face="bold",size=15), 
      legend.text=element_text(face="bold",size=15))
```


    
![png](output_38_0.png)
    



```R
num_cluster = 1

subdata = as.data.frame(UMAP_DIFF_correlation_50_2_1_1[which(UMAP_DIFF_correlation_50_2_1_1$cluster == num_cluster),])

AB_DIFF = c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10")

AB_DIFF = as.data.frame(AB_DIFF)
colnames(AB_DIFF) = "AB_DIFF"
AB_DIFF$score = 0

for(i in 1:nrow(AB_DIFF)){
AB_DIFF[i,2] = nrow(subdata[which(subdata$DIFF_Q == AB_DIFF[i,1]),])
    }

Q10_data = cbind(num_cluster,round(sum(AB_DIFF[10,2])/sum(AB_DIFF[,2])*100,2))
Q10_data = as.data.frame(Q10_data)
colnames(Q10_data) = c("cluster","Q10_ratio")


for(s in 2:22){
num_cluster = s

subdata = as.data.frame(UMAP_DIFF_correlation_50_2_1_1[which(UMAP_DIFF_correlation_50_2_1_1$cluster == num_cluster),])

AB_DIFF = c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10")

AB_DIFF = as.data.frame(AB_DIFF)
colnames(AB_DIFF) = "AB_DIFF"
AB_DIFF$score = 0

for(i in 1:nrow(AB_DIFF)){
AB_DIFF[i,2] = nrow(subdata[which(subdata$DIFF_Q == AB_DIFF[i,1]),])
    }
Q10_data_frag = cbind(num_cluster,round(sum(AB_DIFF[10,2])/sum(AB_DIFF[,2])*100,2))
Q10_data_frag = as.data.frame(Q10_data_frag)
colnames(Q10_data_frag) = c("cluster","Q10_ratio")

    Q10_data = rbind(Q10_data,Q10_data_frag)

    }

Q10_data[,1] = as.factor(Q10_data[,1])
Q10_data$label =  "Q10"
Q10_data = inner_join(DIFF_clustersorted,Q10_data,by="cluster")

Q10_data$sort = c(1:22)
Q10_data[,4] = as.factor(Q10_data[,4])

head(Q10_data)
```


<table class="dataframe">
<caption>A data.frame: 6 × 4</caption>
<thead>
	<tr><th></th><th scope=col>cluster</th><th scope=col>Q10_ratio</th><th scope=col>label</th><th scope=col>sort</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>14</td><td>37.99</td><td>Q10</td><td>1</td></tr>
	<tr><th scope=row>2</th><td>3 </td><td>44.56</td><td>Q10</td><td>2</td></tr>
	<tr><th scope=row>3</th><td>7 </td><td> 0.00</td><td>Q10</td><td>3</td></tr>
	<tr><th scope=row>4</th><td>16</td><td>40.30</td><td>Q10</td><td>4</td></tr>
	<tr><th scope=row>5</th><td>17</td><td> 0.00</td><td>Q10</td><td>5</td></tr>
	<tr><th scope=row>6</th><td>8 </td><td> 1.17</td><td>Q10</td><td>6</td></tr>
</tbody>
</table>




```R
tail(Q10_data,10)
```


<table class="dataframe">
<caption>A data.frame: 10 × 4</caption>
<thead>
	<tr><th></th><th scope=col>cluster</th><th scope=col>Q_ratio</th><th scope=col>label</th><th scope=col>sort</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>13</th><td>6 </td><td> 2.89</td><td>Q10</td><td>13</td></tr>
	<tr><th scope=row>14</th><td>18</td><td> 0.28</td><td>Q10</td><td>14</td></tr>
	<tr><th scope=row>15</th><td>5 </td><td> 3.75</td><td>Q10</td><td>15</td></tr>
	<tr><th scope=row>16</th><td>22</td><td> 0.06</td><td>Q10</td><td>16</td></tr>
	<tr><th scope=row>17</th><td>1 </td><td> 3.40</td><td>Q10</td><td>17</td></tr>
	<tr><th scope=row>18</th><td>9 </td><td> 2.90</td><td>Q10</td><td>18</td></tr>
	<tr><th scope=row>19</th><td>13</td><td> 6.13</td><td>Q10</td><td>19</td></tr>
	<tr><th scope=row>20</th><td>10</td><td> 0.00</td><td>Q10</td><td>20</td></tr>
	<tr><th scope=row>21</th><td>4 </td><td>15.14</td><td>Q10</td><td>21</td></tr>
	<tr><th scope=row>22</th><td>2 </td><td> 0.00</td><td>Q10</td><td>22</td></tr>
</tbody>
</table>




```R
options(repr.plot.width = 5, repr.plot.height = 3, repr.plot.res = 1000, repr.plot.pointsize = 40)
ggplot(Q10_data,aes(sort, y = Q10_ratio, group = "Q10")) +
geom_line(color="darkred", size=0.5, alpha = 1) +
geom_point(color = "black", size = 2) + 
theme_classic(base_size = 10) +
ggtitle(NULL)+
xlab(NULL) +
ylab(NULL) +
theme(plot.title=element_text(face="bold",hjust=0.5,size=25,color = "darkblue"),
      axis.text.x=element_text(face="bold",size=0),
      axis.text.y=element_text(face="bold",size=10),
      axis.title.x = element_text(face="bold",size = 20),
      axis.title.y = element_text(face="bold",size = 15,color = "darkblue"),
      legend.title=element_text(face="bold",size=15), 
      legend.text=element_text(face="bold",size=15))
```


    
![png](output_41_0.png)
    



```R
num_cluster = 1

subdata = as.data.frame(UMAP_DIFF_correlation_50_2_1_1[which(UMAP_DIFF_correlation_50_2_1_1$cluster == num_cluster),])

AB_DIFF = c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10")

AB_DIFF = as.data.frame(AB_DIFF)
colnames(AB_DIFF) = "AB_DIFF"
AB_DIFF$score = 0

for(i in 1:nrow(AB_DIFF)){
AB_DIFF[i,2] = nrow(subdata[which(subdata$DIFF_Q == AB_DIFF[i,1]),])
    }

Q1_data = cbind(num_cluster,round(sum(AB_DIFF[1,2])/sum(AB_DIFF[,2])*100,2))
Q1_data = as.data.frame(Q1_data)
colnames(Q1_data) = c("cluster","Q1_ratio")


for(s in 2:22){
num_cluster = s

subdata = as.data.frame(UMAP_DIFF_correlation_50_2_1_1[which(UMAP_DIFF_correlation_50_2_1_1$cluster == num_cluster),])

AB_DIFF = c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10")

AB_DIFF = as.data.frame(AB_DIFF)
colnames(AB_DIFF) = "AB_DIFF"
AB_DIFF$score = 0

for(i in 1:nrow(AB_DIFF)){
AB_DIFF[i,2] = nrow(subdata[which(subdata$DIFF_Q == AB_DIFF[i,1]),])
    }
Q1_data_frag = cbind(num_cluster,round(sum(AB_DIFF[1,2])/sum(AB_DIFF[,2])*100,2))
Q1_data_frag = as.data.frame(Q1_data_frag)
colnames(Q1_data_frag) = c("cluster","Q1_ratio")

    Q1_data = rbind(Q1_data,Q1_data_frag)

    }

Q1_data[,1] = as.factor(Q1_data[,1])
Q1_data$label =  "Q1"
Q1_data = inner_join(DIFF_clustersorted,Q1_data,by="cluster")

Q1_data$sort = c(1:22)
Q1_data[,4] = as.factor(Q1_data[,4])

head(Q1_data)
```


<table class="dataframe">
<caption>A data.frame: 6 × 4</caption>
<thead>
	<tr><th></th><th scope=col>cluster</th><th scope=col>Q1_ratio</th><th scope=col>label</th><th scope=col>sort</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>14</td><td> 0.00</td><td>Q1</td><td>1</td></tr>
	<tr><th scope=row>2</th><td>3 </td><td>12.15</td><td>Q1</td><td>2</td></tr>
	<tr><th scope=row>3</th><td>7 </td><td> 0.00</td><td>Q1</td><td>3</td></tr>
	<tr><th scope=row>4</th><td>16</td><td> 5.72</td><td>Q1</td><td>4</td></tr>
	<tr><th scope=row>5</th><td>17</td><td> 0.00</td><td>Q1</td><td>5</td></tr>
	<tr><th scope=row>6</th><td>8 </td><td> 0.00</td><td>Q1</td><td>6</td></tr>
</tbody>
</table>




```R
options(repr.plot.width = 5, repr.plot.height = 3, repr.plot.res = 1000, repr.plot.pointsize = 40)
ggplot(Q1_data,aes(sort, y = Q1_ratio, group = "Q1")) +
geom_line(color="darkblue", size=0.5, alpha = 1) +
geom_point(color = "black", size = 2) + 
theme_classic(base_size = 10) +
ggtitle(NULL)+
xlab(NULL) +
ylab(NULL) +
theme(plot.title=element_text(face="bold",hjust=0.5,size=25,color = "darkblue"),
      axis.text.x=element_text(face="bold",size=0),
      axis.text.y=element_text(face="bold",size=10),
      axis.title.x = element_text(face="bold",size = 20),
      axis.title.y = element_text(face="bold",size = 15,color = "darkblue"),
      legend.title=element_text(face="bold",size=15), 
      legend.text=element_text(face="bold",size=15))
```


    
![png](output_43_0.png)
    



```R
colnames(Q1_data) = c("cluster","Q_ratio","label","sort")
colnames(Q10_data) = c("cluster","Q_ratio","label","sort")
QQ_data = rbind(Q1_data,Q10_data)

options(repr.plot.width = 5, repr.plot.height = 3, repr.plot.res = 1000, repr.plot.pointsize = 40)
ggplot(QQ_data,aes(sort, y = Q_ratio,)) +
geom_line(data=subset(QQ_data,label=="Q1"), color="darkblue", size=0.7, alpha = 1, group = "Q1") +
geom_line(data=subset(QQ_data,label=="Q10"), color="darkred", size=0.7, alpha = 1, group = "Q10") +
geom_point(color = "black", size = 2) + 
theme_classic(base_size = 10) +
ggtitle(NULL)+
xlab(NULL) +
ylab(NULL) +
theme(plot.title=element_text(face="bold",hjust=0.5,size=25,color = "darkblue"),
      axis.text.x=element_text(face="bold",size=0),
      axis.text.y=element_text(face="bold",size=10),
      axis.title.x = element_text(face="bold",size = 20),
      axis.title.y = element_text(face="bold",size = 15,color = "darkblue"),
      legend.title=element_text(face="bold",size=15), 
      legend.text=element_text(face="bold",size=15))
```


    
![png](output_44_0.png)
    



```R

```


<table class="dataframe">
<caption>A data.frame: 6 × 4</caption>
<thead>
	<tr><th></th><th scope=col>cluster</th><th scope=col>Q1_ratio</th><th scope=col>label</th><th scope=col>sort</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>14</td><td> 0.00</td><td>Q1</td><td>1</td></tr>
	<tr><th scope=row>2</th><td>3 </td><td>12.15</td><td>Q1</td><td>2</td></tr>
	<tr><th scope=row>3</th><td>7 </td><td> 0.00</td><td>Q1</td><td>3</td></tr>
	<tr><th scope=row>4</th><td>16</td><td> 5.72</td><td>Q1</td><td>4</td></tr>
	<tr><th scope=row>5</th><td>17</td><td> 0.00</td><td>Q1</td><td>5</td></tr>
	<tr><th scope=row>6</th><td>8 </td><td> 0.00</td><td>Q1</td><td>6</td></tr>
</tbody>
</table>




```R

```


```R

```


```R
# 14 / 13
num_cluster = 9

subdata = as.data.frame(UMAP_DIFF_correlation_50_2_1_1[which(UMAP_DIFF_correlation_50_2_1_1$cluster == num_cluster),])

Ets1_Q = c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10")

Ets1_Q = as.data.frame(Ets1_Q)
colnames(Ets1_Q) = "Ets1_Q"
Ets1_Q$score = 0

for(i in 1:nrow(Ets1_Q)){
Ets1_Q[i,2] = nrow(subdata[which(subdata$Ets1_Q == Ets1_Q[i,1]),])
    }

#Ets1_Q = Ets1_Q[order(-Ets1_Q$score),]

subdata = as.data.frame(UMAP_DIFF_correlation_50_2_1_1[which(UMAP_DIFF_correlation_50_2_1_1$cluster == num_cluster),])

p300_Q = c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10")

p300_Q = as.data.frame(p300_Q)
colnames(p300_Q) = "p300_Q"
p300_Q$score = 0

for(i in 1:nrow(p300_Q)){
p300_Q[i,2] = nrow(subdata[which(subdata$p300_Q == p300_Q[i,1]),])
    }

#p300_Q = p300_Q[order(-p300_Q$score),]

ddata= cbind(Ets1_Q,p300_Q)

ddata$label = c("A","B","C","D","E","F","G","H","I","J")

ddata
```


<table class="dataframe">
<caption>A data.frame: 10 × 5</caption>
<thead>
	<tr><th scope=col>Ets1_Q</th><th scope=col>score</th><th scope=col>p300_Q</th><th scope=col>score</th><th scope=col>label</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Q1 </td><td>264</td><td>Q1 </td><td>246</td><td>A</td></tr>
	<tr><td>Q2 </td><td>128</td><td>Q2 </td><td>152</td><td>B</td></tr>
	<tr><td>Q3 </td><td>116</td><td>Q3 </td><td>120</td><td>C</td></tr>
	<tr><td>Q4 </td><td> 91</td><td>Q4 </td><td>108</td><td>D</td></tr>
	<tr><td>Q5 </td><td>104</td><td>Q5 </td><td>100</td><td>E</td></tr>
	<tr><td>Q6 </td><td>101</td><td>Q6 </td><td>100</td><td>F</td></tr>
	<tr><td>Q7 </td><td> 91</td><td>Q7 </td><td> 80</td><td>G</td></tr>
	<tr><td>Q8 </td><td> 82</td><td>Q8 </td><td> 71</td><td>H</td></tr>
	<tr><td>Q9 </td><td>104</td><td>Q9 </td><td> 95</td><td>I</td></tr>
	<tr><td>Q10</td><td> 90</td><td>Q10</td><td> 99</td><td>J</td></tr>
</tbody>
</table>




```R
options(repr.plot.width = 8, repr.plot.height = 8, repr.plot.res = 1000, repr.plot.pointsize = 10)

data <- data.frame(
    category=ddata[,5],
    count=ddata[,2])
 
data$fraction <- data$count / sum(data$count)
data$ymax <- cumsum(data$fraction)
data$ymin <- c(0, head(data$ymax, n=-1))
data$labelPosition <- (data$ymax + data$ymin) / 2


data$label <- paste0(data$category, "\n value: ", data$count)
nb.cols <- 19
mycolors <- colorRampPalette(brewer.pal(8, "Pastel2"))(nb.cols)


ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
geom_rect() +
scale_color_manual(values=c("grey40"))+
scale_fill_manual(values=c("#FFF8E1","#FFECB3","#FFE082","#FFD54F","#FFCA28",
                           "#FFC107","#FFB300","#FFA000","#FF8F00","#FF6F00"))+
coord_polar(theta="y") +
  coord_polar(theta="y") +
  xlim(c(2.5, 4)) +
  theme_void() +
scale_colour_hue(l = 20, c = 100)+
  theme(legend.position = "none")

data <- data.frame(
    category=ddata[,5],
    count=ddata[,4])
 
data$fraction <- data$count / sum(data$count)
data$ymax <- cumsum(data$fraction)
data$ymin <- c(0, head(data$ymax, n=-1))
data$labelPosition <- (data$ymax + data$ymin) / 2


data$label <- paste0(data$category, "\n value: ", data$count)
nb.cols <- 19
mycolors <- colorRampPalette(brewer.pal(8, "Pastel2"))(nb.cols)


ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
geom_rect() +
scale_color_manual(values=c("grey40"))+
scale_fill_manual(values=c("#E0F2F1","#B2DFDB","#80CBC4","#4DB6AC","#26A69A",
                           "#009688","#00897B","#00796B","#00695C","#004D40"))+
coord_polar(theta="y") +
  coord_polar(theta="y") +
  xlim(c(2.5, 4)) +
  theme_void() +
scale_colour_hue(l = 20, c = 100)+
  theme(legend.position = "none")



```

    [1m[22mCoordinate system already present. Adding new coordinate system, which will replace the existing one.
    [1m[22mScale for [32mcolour[39m is already present.
    Adding another scale for [32mcolour[39m, which will replace the existing scale.
    [1m[22mCoordinate system already present. Adding new coordinate system, which will replace the existing one.
    [1m[22mScale for [32mcolour[39m is already present.
    Adding another scale for [32mcolour[39m, which will replace the existing scale.



    
![png](output_49_1.png)
    



    
![png](output_49_2.png)
    



```R
options(repr.plot.width = 8, repr.plot.height = 8, repr.plot.res = 1000, repr.plot.pointsize = 10)

data <- data.frame(
    category=ddata[,5],
    count=ddata[,2])
 
data$fraction <- data$count / sum(data$count)
data$ymax <- cumsum(data$fraction)
data$ymin <- c(0, head(data$ymax, n=-1))
data$labelPosition <- (data$ymax + data$ymin) / 2


data$label <- paste0(data$category, "\n value: ", data$count)
nb.cols <- 19
mycolors <- colorRampPalette(brewer.pal(8, "Pastel2"))(nb.cols)


ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
geom_rect() +
scale_color_manual(values=c("grey40"))+
scale_fill_manual(values=c("#FFF8E1","#FFECB3","#FFE082","#FFD54F","#FFCA28",
                           "#FFC107","#FFB300","#FFA000","#FF8F00","#FF6F00"))+
coord_polar(theta="y") +
  coord_polar(theta="y") +
  xlim(c(2.5, 4)) +
  theme_void() +
scale_colour_hue(l = 20, c = 100)+
  theme(legend.position = "none")

data <- data.frame(
    category=ddata[,5],
    count=ddata[,4])
 
data$fraction <- data$count / sum(data$count)
data$ymax <- cumsum(data$fraction)
data$ymin <- c(0, head(data$ymax, n=-1))
data$labelPosition <- (data$ymax + data$ymin) / 2


data$label <- paste0(data$category, "\n value: ", data$count)
nb.cols <- 19
mycolors <- colorRampPalette(brewer.pal(8, "Pastel2"))(nb.cols)


ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
geom_rect() +
scale_color_manual(values=c("grey40"))+
scale_fill_manual(values=c("#E0F2F1","#B2DFDB","#80CBC4","#4DB6AC","#26A69A",
                           "#009688","#00897B","#00796B","#00695C","#004D40"))+
coord_polar(theta="y") +
  coord_polar(theta="y") +
  xlim(c(2.5, 4)) +
  theme_void() +
scale_colour_hue(l = 20, c = 100)+
  theme(legend.position = "none")



```

    [1m[22mCoordinate system already present. Adding new coordinate system, which will replace the existing one.
    [1m[22mScale for [32mcolour[39m is already present.
    Adding another scale for [32mcolour[39m, which will replace the existing scale.
    [1m[22mCoordinate system already present. Adding new coordinate system, which will replace the existing one.
    [1m[22mScale for [32mcolour[39m is already present.
    Adding another scale for [32mcolour[39m, which will replace the existing scale.



    
![png](output_50_1.png)
    



    
![png](output_50_2.png)
    



```R

```


```R

```


```R

```


```R

```


```R


i = 1

ABsignal_cluster = cbind(i,median(UMAP_DIFF_correlation_50_2_1_1[which(UMAP_DIFF_correlation_50_2_1_1$cluster == i),8]))
ABsignal_cluster = as.data.frame(ABsignal_cluster)
colnames(ABsignal_cluster) = c("cluster","ABsignal")

for(i in 2:22){
    ABsignal_cluster_frag = cbind(i,median(UMAP_DIFF_correlation_50_2_1_1[which(UMAP_DIFF_correlation_50_2_1_1$cluster == i),8]))
    ABsignal_cluster_frag = as.data.frame(ABsignal_cluster_frag)
    colnames(ABsignal_cluster_frag) = c("cluster","ABsignal")

    ABsignal_cluster = rbind(ABsignal_cluster,ABsignal_cluster_frag)
    }

ABsignal_cluster = ABsignal_cluster[order(-ABsignal_cluster$ABsignal),]
rownames(ABsignal_cluster) = NULL
ABsignal_cluster[,1] = as.factor(ABsignal_cluster[,1])

UMAP_DIFF_correlation_50_2_1_1_forPlot_DIFF = UMAP_DIFF_correlation_50_2_1_1[,c(14,8)]

UMAP_DIFF_correlation_50_2_1_1_forPlot_Ets1 = UMAP_DIFF_correlation_50_2_1_1[,c(14,10)]

UMAP_DIFF_correlation_50_2_1_1_forPlot_Ets1[,2] = log2(UMAP_DIFF_correlation_50_2_1_1_forPlot_Ets1[,2]+1)

head(UMAP_DIFF_correlation_50_2_1_1)
```


```R

```


```R

```


<table class="dataframe">
<caption>A data.frame: 6 × 16</caption>
<thead>
	<tr><th></th><th scope=col>V1</th><th scope=col>V2</th><th scope=col>V3</th><th scope=col>V4</th><th scope=col>V5</th><th scope=col>V6</th><th scope=col>V7</th><th scope=col>V8</th><th scope=col>V9</th><th scope=col>V10</th><th scope=col>V11</th><th scope=col>V12</th><th scope=col>V13</th><th scope=col>V14</th><th scope=col>V15</th><th scope=col>V16</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>chr1</td><td>3000001</td><td>3100000</td><td> 0.0005551651</td><td>-0.004033565</td><td>-0.004588731</td><td>Q4</td><td> 52</td><td>Q1 </td><td>68</td><td>Q6 </td><td>1</td><td>chr1</td><td>3072253</td><td>3074252</td><td>4933401J01Rik</td></tr>
	<tr><th scope=row>2</th><td>chr1</td><td>3100001</td><td>3200000</td><td> 0.0005551651</td><td>-0.003834151</td><td>-0.004389316</td><td>Q4</td><td> 75</td><td>Q4 </td><td>77</td><td>Q8 </td><td>1</td><td>chr1</td><td>3101016</td><td>3103015</td><td>Gm26206      </td></tr>
	<tr><th scope=row>3</th><td>chr1</td><td>3200001</td><td>3300000</td><td> 0.0005551651</td><td>-0.003926483</td><td>-0.004481648</td><td>Q4</td><td> 90</td><td>Q8 </td><td>62</td><td>Q4 </td><td>1</td><td>chr1</td><td>3251757</td><td>3253756</td><td>Gm18956      </td></tr>
	<tr><th scope=row>4</th><td>chr1</td><td>3300001</td><td>3400000</td><td>-0.0001989466</td><td>-0.004390436</td><td>-0.004191489</td><td>Q4</td><td>102</td><td>Q10</td><td>82</td><td>Q9 </td><td>1</td><td>chr1</td><td>3367550</td><td>3369549</td><td>Gm37180      </td></tr>
	<tr><th scope=row>5</th><td>chr1</td><td>3300001</td><td>3400000</td><td>-0.0001989466</td><td>-0.004390436</td><td>-0.004191489</td><td>Q4</td><td>102</td><td>Q10</td><td>82</td><td>Q9 </td><td>1</td><td>chr1</td><td>3376789</td><td>3378788</td><td>Gm37363      </td></tr>
	<tr><th scope=row>6</th><td>chr1</td><td>3400001</td><td>3500000</td><td> 0.0020771037</td><td>-0.005041151</td><td>-0.007118255</td><td>Q4</td><td>105</td><td>Q10</td><td>94</td><td>Q10</td><td>1</td><td>chr1</td><td>3465587</td><td>3467586</td><td>Gm1992       </td></tr>
</tbody>
</table>




```R

```


```R

```


```R

```


```R

```


```R

```


```R

```


```R

```


```R
library(GenomicFeatures)
library(rtracklayer)


gtf <- import("/data3/psg/gencode_mm10_annotation.gtf")
gtf_genes <- gtf[gtf$type == "gene"]
gene_info <- mcols(gtf_genes)[, c("gene_id", "gene_name")]
txdb <- makeTxDbFromGFF("/data3/psg/gencode_mm10_annotation.gtf")
promoters_gr <- promoters(genes(txdb), upstream=1000, downstream=1000)
promoters_df <- as.data.frame(promoters_gr)

promoters_df$gene_id <- rownames(promoters_df)  
merged <- merge(promoters_df, gene_info, by.x="gene_id", by.y="gene_id", all.x=TRUE)

saf <- data.frame(
  GeneID = merged$gene_name,
  Chr = merged$seqnames,
  Start = merged$start,
  End = merged$end,
  Strand = merged$strand
)

saf = saf[which(saf$Start > 0),]
saf = saf[which(saf$End > 0),]

saf =  saf[,c(2,3,4,1)]

saf=saf[!duplicated(saf[,4]),]
saf = saf[grepl("chr",saf[,1]), ]


write.table(saf, "/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/DATAanalysis/gencode_mm10_annotation_promoters.bed",
            sep="\t", quote=FALSE, row.names=FALSE)
```

    Import genomic features from the file as a GRanges object ... 
    OK
    
    Prepare the 'metadata' data frame ... 
    OK
    
    Make the TxDb object ... 
    Warning message in .get_cds_IDX(mcols0$type, mcols0$phase):
    “The "phase" metadata column contains non-NA values for features of type
      stop_codon. This information was ignored.”
    OK
    



```R
set.seed(7777)

UMAP_DIFF_correlation_50_2_1_1 = as.data.frame(vroom("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_50_2_1_1.txt", delim = "\t",show_col_types = FALSE))
UMAP_DIFF_correlation_50_2_1_1_fit = UMAP_DIFF_correlation_50_2_1_1[,1:2]
UMAP_DIFF_correlation_50_2_1_1_label = paste(UMAP_DIFF_correlation_50_2_1_1[,3],UMAP_DIFF_correlation_50_2_1_1[,4],UMAP_DIFF_correlation_50_2_1_1[,5])
UMAP_DIFF_correlation_50_2_1_1_label = as.data.frame(UMAP_DIFF_correlation_50_2_1_1_label)
colnames(UMAP_DIFF_correlation_50_2_1_1_label) = "label"

result <- kmeans(UMAP_DIFF_correlation_50_2_1_1_fit, centers = 22, iter.max = 10000, algorithm = "MacQueen")

UMAP_DIFF_correlation_50_2_1_1 = cbind(UMAP_DIFF_correlation_50_2_1_1,as.data.frame(result$cluster))
colnames(UMAP_DIFF_correlation_50_2_1_1)[ncol(UMAP_DIFF_correlation_50_2_1_1)] = "cluster"

UMAP_DIFF_correlation_50_2_1_1_fit = UMAP_DIFF_correlation_50_2_1_1[,3:14]

write.table(UMAP_DIFF_correlation_50_2_1_1_fit, 
            file = "/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/UMAP_DIFF_correlation_50_2_1_1_fit.txt" ,col.names=F, row.names=F, quote=F, sep="\t")

```


```R
# C14 Igf1
```


```R
UMAP_DIFF_correlation_50_2_1_1_fittoPROMOTER = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/DATAanalysis/UMAP_DIFF_correlation_50_2_1_1_fittoPROMOTER.txt",
                                                          sep = "\t", header=F,fill = TRUE)

UMAP_DIFF_correlation_50_2_1_1_fittoPROMOTER[which(UMAP_DIFF_correlation_50_2_1_1_fittoPROMOTER[,16] == "Tcf21"),]

C14 = UMAP_DIFF_correlation_50_2_1_1_fittoPROMOTER[which(UMAP_DIFF_correlation_50_2_1_1_fittoPROMOTER[,12] == 14),]
C14 = C14[grepl("Q10|Q9|Q8",C14[,7]), ]
C14 = C14[grepl("Q10|Q9|Q8",C14[,9]), ]
C14 = C14[grepl("Q10|Q9|Q8",C14[,11]), ]
C14 = C14[!grepl("Gm",C14[,16]), ]
C14 = C14[!grepl("Rik",C14[,16]), ]

C14[,16]
```


<table class="dataframe">
<caption>A data.frame: 1 × 16</caption>
<thead>
	<tr><th></th><th scope=col>V1</th><th scope=col>V2</th><th scope=col>V3</th><th scope=col>V4</th><th scope=col>V5</th><th scope=col>V6</th><th scope=col>V7</th><th scope=col>V8</th><th scope=col>V9</th><th scope=col>V10</th><th scope=col>V11</th><th scope=col>V12</th><th scope=col>V13</th><th scope=col>V14</th><th scope=col>V15</th><th scope=col>V16</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>3677</th><td>chr10</td><td>22800001</td><td>22900000</td><td>-0.0176519</td><td>0.02899706</td><td>0.04664896</td><td>Q10</td><td>89</td><td>Q8</td><td>77</td><td>Q8</td><td>14</td><td>chr10</td><td>22819129</td><td>22821128</td><td>Tcf21</td></tr>
</tbody>
</table>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Mthfd1l'</li><li>'Stxbp5'</li><li>'Grm1'</li><li>'Perp'</li><li>'Olig3'</li><li>'Il20ra'</li><li>'Tbpl1'</li><li>'Mir7663'</li><li>'Tcf21'</li><li>'Cenpw'</li><li>'Tpd52l1'</li><li>'Nt5dc1'</li><li>'Frk'</li><li>'Hs3st5'</li><li>'Prdm1'</li><li>'Grik2'</li><li>'Asf1a'</li><li>'Mcm9'</li><li>'Man1a'</li><li>'Msl3l2'</li><li>'Hsf2'</li><li>'AW822073'</li><li>'Duxf3'</li><li>'Gcc2'</li><li>'Lrrtm3'</li><li>'Olfr1354'</li><li>'Olfr8'</li><li>'Olfr1353'</li><li>'Olfr1352'</li><li>'Zfp781'</li><li>'Timp3'</li><li>'Stab2'</li><li>'Ascl1'</li><li>'Pah'</li><li>'Igf1'</li><li>'Igf1os'</li><li>'Pmch'</li><li>'Chpt1'</li><li>'Ano4'</li><li>'Mir135a-2'</li><li>'Mir1251'</li><li>'Rmst'</li><li>'Cfap54'</li><li>'Rassf9'</li><li>'Alx1'</li><li>'Tmtc2'</li><li>'Acss3'</li><li>'Lin7a'</li><li>'Myf5'</li><li>'Myf6'</li><li>'Ptprq'</li><li>'Nav3'</li><li>'Nav3'</li><li>'E2f7'</li><li>'Tph2'</li><li>'Tspan8'</li><li>'Grip1os2'</li><li>'Mon2'</li><li>'Usp15'</li><li>'Lrig3'</li><li>'Olfr786'</li><li>'Olfr787'</li><li>'Olfr788'</li><li>'Olfr789'</li><li>'Olfr790'</li><li>'Olfr827'</li><li>'Neurod4'</li><li>'Eif3s6-ps1'</li><li>'Vwc2'</li><li>'Ikzf1'</li><li>'Bcl11a'</li><li>'Anp32-ps'</li><li>'Tenm2'</li><li>'Hmgb1-ps1'</li><li>'Gabrg2'</li><li>'Gabra1'</li><li>'Gabrb2'</li></ol>




```R
UMAP_DIFF_correlation_50_2_1_1_fittoPROMOTER = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/DATAanalysis/UMAP_DIFF_correlation_50_2_1_1_fittoPROMOTER.txt",
                                                          sep = "\t", header=F,fill = TRUE)

UMAP_DIFF_correlation_50_2_1_1_fittoPROMOTER[which(UMAP_DIFF_correlation_50_2_1_1_fittoPROMOTER[,16] == "Sod2"),]

C14 = UMAP_DIFF_correlation_50_2_1_1_fittoPROMOTER[which(UMAP_DIFF_correlation_50_2_1_1_fittoPROMOTER[,12] == 16),]
C14 = C14[grepl("Q10|Q9|Q8",C14[,7]), ]
C14 = C14[grepl("Q10|Q9|Q8",C14[,9]), ]
C14 = C14[grepl("Q10|Q9|Q8",C14[,11]), ]
C14 = C14[!grepl("Gm",C14[,16]), ]
C14 = C14[!grepl("Rik",C14[,16]), ]
```


<table class="dataframe">
<caption>A data.frame: 1 × 16</caption>
<thead>
	<tr><th></th><th scope=col>V1</th><th scope=col>V2</th><th scope=col>V3</th><th scope=col>V4</th><th scope=col>V5</th><th scope=col>V6</th><th scope=col>V7</th><th scope=col>V8</th><th scope=col>V9</th><th scope=col>V10</th><th scope=col>V11</th><th scope=col>V12</th><th scope=col>V13</th><th scope=col>V14</th><th scope=col>V15</th><th scope=col>V16</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>15800</th><td>chr17</td><td>13000001</td><td>13100000</td><td>-0.02443585</td><td>0.01211657</td><td>0.03655242</td><td>Q10</td><td>109</td><td>Q10</td><td>119</td><td>Q10</td><td>16</td><td>chr17</td><td>13006839</td><td>13008838</td><td>Sod2</td></tr>
</tbody>
</table>




```R
UMAP_DIFF_correlation_50_2_1_1_fittoPROMOTER = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/DATAanalysis/UMAP_DIFF_correlation_50_2_1_1_fittoPROMOTER.txt",
                                                          sep = "\t", header=F,fill = TRUE)

UMAP_DIFF_correlation_50_2_1_1_fittoPROMOTER[which(UMAP_DIFF_correlation_50_2_1_1_fittoPROMOTER[,16] == "Rrm2"),]

C14 = UMAP_DIFF_correlation_50_2_1_1_fittoPROMOTER[which(UMAP_DIFF_correlation_50_2_1_1_fittoPROMOTER[,12] == 13),]
C14 = C14[grepl("Q1$|Q2|Q3",C14[,7]), ]
C14 = C14[grepl("Q1$|Q2|Q3",C14[,9]), ]
C14 = C14[grepl("Q1$|Q2|Q3",C14[,11]), ]
C14 = C14[!grepl("Gm",C14[,16]), ]
C14 = C14[!grepl("Rik",C14[,16]), ]
```


<table class="dataframe">
<caption>A data.frame: 1 × 16</caption>
<thead>
	<tr><th></th><th scope=col>V1</th><th scope=col>V2</th><th scope=col>V3</th><th scope=col>V4</th><th scope=col>V5</th><th scope=col>V6</th><th scope=col>V7</th><th scope=col>V8</th><th scope=col>V9</th><th scope=col>V10</th><th scope=col>V11</th><th scope=col>V12</th><th scope=col>V13</th><th scope=col>V14</th><th scope=col>V15</th><th scope=col>V16</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>8394</th><td>chr12</td><td>24700001</td><td>24800000</td><td>0.01308099</td><td>-0.01423885</td><td>-0.02731983</td><td>Q2</td><td>62</td><td>Q2</td><td>44</td><td>Q1</td><td>13</td><td>chr12</td><td>24707241</td><td>24709240</td><td>Rrm2</td></tr>
</tbody>
</table>




```R
UMAP_DIFF_correlation_50_2_1_1_fittoPROMOTER = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG5/UMAP/DATAanalysis/UMAP_DIFF_correlation_50_2_1_1_fittoPROMOTER.txt",
                                                          sep = "\t", header=F,fill = TRUE)

UMAP_DIFF_correlation_50_2_1_1_fittoPROMOTER[which(UMAP_DIFF_correlation_50_2_1_1_fittoPROMOTER[,16] == "Usp1"),]

C14 = UMAP_DIFF_correlation_50_2_1_1_fittoPROMOTER[which(UMAP_DIFF_correlation_50_2_1_1_fittoPROMOTER[,12] == 9),]
C14 = C14[grepl("Q1$|Q2|Q3",C14[,7]), ]
C14 = C14[grepl("Q1$|Q2|Q3",C14[,9]), ]
C14 = C14[grepl("Q1$|Q2|Q3",C14[,11]), ]
C14 = C14[!grepl("Gm",C14[,16]), ]
C14 = C14[!grepl("Rik",C14[,16]), ]
```


<table class="dataframe">
<caption>A data.frame: 1 × 16</caption>
<thead>
	<tr><th></th><th scope=col>V1</th><th scope=col>V2</th><th scope=col>V3</th><th scope=col>V4</th><th scope=col>V5</th><th scope=col>V6</th><th scope=col>V7</th><th scope=col>V8</th><th scope=col>V9</th><th scope=col>V10</th><th scope=col>V11</th><th scope=col>V12</th><th scope=col>V13</th><th scope=col>V14</th><th scope=col>V15</th><th scope=col>V16</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>27543</th><td>chr4</td><td>98900001</td><td>99000000</td><td>0.03374925</td><td>-0.004789652</td><td>-0.0385389</td><td>Q2</td><td>65</td><td>Q2</td><td>50</td><td>Q2</td><td>9</td><td>chr4</td><td>98922810</td><td>98924809</td><td>Usp1</td></tr>
</tbody>
</table>




```R

```


```R

```


```R

```


```R

```


```R

```


```R

```


```R

```


```R

```


```R

```
