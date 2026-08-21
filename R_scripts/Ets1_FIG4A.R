ATAC_Full_data = read.table("/data3/psg/Ets1/ProcessedData/Diff_Access_TF_Full_data.txt",sep = "\t", header=T,fill = TRUE)
TFBS_Full_data = read.table("/data3/psg/Ets1/ProcessedData/Diff_TFBS_TF_Full_data.txt",sep = "\t", header=T,fill = TRUE)
ATAC_DESeq2_Norm = read.table("/data3/psg/Ets1/ProcessedData/Access_TF_count_data.txt",sep = "\t", header=T,fill = TRUE)
TFBS_DESeq2_Norm = read.table("/data3/psg/Ets1/ProcessedData/TFBS_TF_count_data.txt",sep = "\t", header=T,fill = TRUE)

ATAC_UP = ATAC_Full_data[which(ATAC_Full_data$log2FoldChange > log2(1.2)),]
ATAC_UP = ATAC_UP[which(ATAC_UP$pvalue < 0.05),]
ATAC_UP = ATAC_UP[order(-ATAC_UP$log2FoldChange),]

ATAC_DOWN = ATAC_Full_data[which(ATAC_Full_data$log2FoldChange < -log2(1.2)),]
ATAC_DOWN = ATAC_DOWN[which(ATAC_DOWN$pvalue < 0.05),]

TFBS_UP = TFBS_Full_data[which(TFBS_Full_data$log2FoldChange > log2(1.2)),]
TFBS_UP = TFBS_UP[which(TFBS_UP$pvalue < 0.05),]
TFBS_UP = TFBS_UP[order(-TFBS_UP$log2FoldChange),]

TFBS_DOWN = TFBS_Full_data[which(TFBS_Full_data$log2FoldChange < -log2(1.2)),]
TFBS_DOWN = TFBS_DOWN[which(TFBS_DOWN$pvalue < 0.05),]
TFBS_DOWN = TFBS_DOWN[order(-TFBS_DOWN$log2FoldChange),]


ATAC_UP_TF = ATAC_UP[,1]
ATAC_UP_TF = as.data.frame(ATAC_UP_TF)
colnames(ATAC_UP_TF) = "gene"

ATAC_DOWN_TF = ATAC_DOWN[,1]
ATAC_DOWN_TF = as.data.frame(ATAC_DOWN_TF)
colnames(ATAC_DOWN_TF) = "gene"

TFBS_UP_TF = TFBS_UP[,1]
TFBS_UP_TF = as.data.frame(TFBS_UP_TF)
colnames(TFBS_UP_TF) = "gene"

TFBS_DOWN_TF = TFBS_DOWN[,1]
TFBS_DOWN_TF = as.data.frame(TFBS_DOWN_TF)
colnames(TFBS_DOWN_TF) = "gene"

ATAC_UP_TFBS_UP_TF = inner_join(ATAC_UP_TF,TFBS_UP_TF,by="gene")
ATAC_UP_TFBS_DOWN_TF = inner_join(ATAC_UP_TF,TFBS_DOWN_TF,by="gene")

ATAC_STEADY_TFBS_UP_TF = anti_join(TFBS_UP_TF,ATAC_UP_TFBS_UP_TF,by="gene")
ATAC_STEADY_TFBS_DOWN_TF = anti_join(TFBS_DOWN_TF,ATAC_UP_TFBS_DOWN_TF,by="gene")

nrow(ATAC_UP_TFBS_UP_TF)+
nrow(ATAC_UP_TFBS_DOWN_TF)+
nrow(ATAC_STEADY_TFBS_UP_TF)+
nrow(ATAC_STEADY_TFBS_DOWN_TF)

nrow(ATAC_UP_TFBS_UP_TF)
nrow(ATAC_UP_TFBS_DOWN_TF)
nrow(ATAC_STEADY_TFBS_UP_TF)
nrow(ATAC_STEADY_TFBS_DOWN_TF)

ATAC_Z_Norm = ATAC_DESeq2_Norm

ATAC_DESeq2_Norm1 = ATAC_DESeq2_Norm[,2]
ATAC_DESeq2_Norm2 = ATAC_DESeq2_Norm[,3]
ATAC_DESeq2_Norm3 = ATAC_DESeq2_Norm[,4]
ATAC_DESeq2_Norm4 = ATAC_DESeq2_Norm[,5]
ATAC_DESeq2_Norm5 = ATAC_DESeq2_Norm[,6]
ATAC_DESeq2_Norm6 = ATAC_DESeq2_Norm[,7]
ATAC_DESeq2_Norm7 = ATAC_DESeq2_Norm[,8]
ATAC_DESeq2_Norm8 = ATAC_DESeq2_Norm[,9]

ATAC_DESeq2_Norm1 = as.data.frame(ATAC_DESeq2_Norm1)
ATAC_DESeq2_Norm2 = as.data.frame(ATAC_DESeq2_Norm2)
ATAC_DESeq2_Norm3 = as.data.frame(ATAC_DESeq2_Norm3)
ATAC_DESeq2_Norm4 = as.data.frame(ATAC_DESeq2_Norm4)
ATAC_DESeq2_Norm5 = as.data.frame(ATAC_DESeq2_Norm5)
ATAC_DESeq2_Norm6 = as.data.frame(ATAC_DESeq2_Norm6)
ATAC_DESeq2_Norm7 = as.data.frame(ATAC_DESeq2_Norm7)
ATAC_DESeq2_Norm8 = as.data.frame(ATAC_DESeq2_Norm8)

colnames(ATAC_DESeq2_Norm1) = "value"
colnames(ATAC_DESeq2_Norm2) = "value"
colnames(ATAC_DESeq2_Norm3) = "value"
colnames(ATAC_DESeq2_Norm4) = "value"
colnames(ATAC_DESeq2_Norm5) = "value"
colnames(ATAC_DESeq2_Norm6) = "value"
colnames(ATAC_DESeq2_Norm7) = "value"
colnames(ATAC_DESeq2_Norm8) = "value"

ATAC_DESeq2_Norm_vec = rbind(ATAC_DESeq2_Norm1,
                             ATAC_DESeq2_Norm2,
                             ATAC_DESeq2_Norm3,
                             ATAC_DESeq2_Norm4,
                             ATAC_DESeq2_Norm5,
                             ATAC_DESeq2_Norm6,
                             ATAC_DESeq2_Norm7,
                             ATAC_DESeq2_Norm8)

ATAC_DESeq2_Norm_vec$norm = (ATAC_DESeq2_Norm_vec[,1] - mean(ATAC_DESeq2_Norm_vec[,1]))/(sd(ATAC_DESeq2_Norm_vec[,1]))

ATAC_Z_Norm[,2] = ATAC_DESeq2_Norm_vec[(1+(nrow(ATAC_Z_Norm)*0)):(nrow(ATAC_Z_Norm)*1),2]
ATAC_Z_Norm[,3] = ATAC_DESeq2_Norm_vec[(1+(nrow(ATAC_Z_Norm)*1)):(nrow(ATAC_Z_Norm)*2),2]
ATAC_Z_Norm[,4] = ATAC_DESeq2_Norm_vec[(1+(nrow(ATAC_Z_Norm)*2)):(nrow(ATAC_Z_Norm)*3),2]
ATAC_Z_Norm[,5] = ATAC_DESeq2_Norm_vec[(1+(nrow(ATAC_Z_Norm)*3)):(nrow(ATAC_Z_Norm)*4),2]
ATAC_Z_Norm[,6] = ATAC_DESeq2_Norm_vec[(1+(nrow(ATAC_Z_Norm)*4)):(nrow(ATAC_Z_Norm)*5),2]
ATAC_Z_Norm[,7] = ATAC_DESeq2_Norm_vec[(1+(nrow(ATAC_Z_Norm)*5)):(nrow(ATAC_Z_Norm)*6),2]
ATAC_Z_Norm[,8] = ATAC_DESeq2_Norm_vec[(1+(nrow(ATAC_Z_Norm)*6)):(nrow(ATAC_Z_Norm)*7),2]
ATAC_Z_Norm[,9] = ATAC_DESeq2_Norm_vec[(1+(nrow(ATAC_Z_Norm)*7)):(nrow(ATAC_Z_Norm)*8),2]

TFBS_Z_Norm = TFBS_DESeq2_Norm

TFBS_DESeq2_Norm1 = TFBS_DESeq2_Norm[,2]
TFBS_DESeq2_Norm2 = TFBS_DESeq2_Norm[,3]
TFBS_DESeq2_Norm3 = TFBS_DESeq2_Norm[,4]
TFBS_DESeq2_Norm4 = TFBS_DESeq2_Norm[,5]
TFBS_DESeq2_Norm5 = TFBS_DESeq2_Norm[,6]
TFBS_DESeq2_Norm6 = TFBS_DESeq2_Norm[,7]
TFBS_DESeq2_Norm7 = TFBS_DESeq2_Norm[,8]
TFBS_DESeq2_Norm8 = TFBS_DESeq2_Norm[,9]

TFBS_DESeq2_Norm1 = as.data.frame(TFBS_DESeq2_Norm1)
TFBS_DESeq2_Norm2 = as.data.frame(TFBS_DESeq2_Norm2)
TFBS_DESeq2_Norm3 = as.data.frame(TFBS_DESeq2_Norm3)
TFBS_DESeq2_Norm4 = as.data.frame(TFBS_DESeq2_Norm4)
TFBS_DESeq2_Norm5 = as.data.frame(TFBS_DESeq2_Norm5)
TFBS_DESeq2_Norm6 = as.data.frame(TFBS_DESeq2_Norm6)
TFBS_DESeq2_Norm7 = as.data.frame(TFBS_DESeq2_Norm7)
TFBS_DESeq2_Norm8 = as.data.frame(TFBS_DESeq2_Norm8)

colnames(TFBS_DESeq2_Norm1) = "value"
colnames(TFBS_DESeq2_Norm2) = "value"
colnames(TFBS_DESeq2_Norm3) = "value"
colnames(TFBS_DESeq2_Norm4) = "value"
colnames(TFBS_DESeq2_Norm5) = "value"
colnames(TFBS_DESeq2_Norm6) = "value"
colnames(TFBS_DESeq2_Norm7) = "value"
colnames(TFBS_DESeq2_Norm8) = "value"

TFBS_DESeq2_Norm_vec = rbind(TFBS_DESeq2_Norm1,
                             TFBS_DESeq2_Norm2,
                             TFBS_DESeq2_Norm3,
                             TFBS_DESeq2_Norm4,
                             TFBS_DESeq2_Norm5,
                             TFBS_DESeq2_Norm6,
                             TFBS_DESeq2_Norm7,
                             TFBS_DESeq2_Norm8)

TFBS_DESeq2_Norm_vec$norm = (TFBS_DESeq2_Norm_vec[,1] - mean(TFBS_DESeq2_Norm_vec[,1]))/(sd(TFBS_DESeq2_Norm_vec[,1]))

TFBS_Z_Norm[,2] = TFBS_DESeq2_Norm_vec[(1+(nrow(TFBS_Z_Norm)*0)):(nrow(TFBS_Z_Norm)*1),2]
TFBS_Z_Norm[,3] = TFBS_DESeq2_Norm_vec[(1+(nrow(TFBS_Z_Norm)*1)):(nrow(TFBS_Z_Norm)*2),2]
TFBS_Z_Norm[,4] = TFBS_DESeq2_Norm_vec[(1+(nrow(TFBS_Z_Norm)*2)):(nrow(TFBS_Z_Norm)*3),2]
TFBS_Z_Norm[,5] = TFBS_DESeq2_Norm_vec[(1+(nrow(TFBS_Z_Norm)*3)):(nrow(TFBS_Z_Norm)*4),2]
TFBS_Z_Norm[,6] = TFBS_DESeq2_Norm_vec[(1+(nrow(TFBS_Z_Norm)*4)):(nrow(TFBS_Z_Norm)*5),2]
TFBS_Z_Norm[,7] = TFBS_DESeq2_Norm_vec[(1+(nrow(TFBS_Z_Norm)*5)):(nrow(TFBS_Z_Norm)*6),2]
TFBS_Z_Norm[,8] = TFBS_DESeq2_Norm_vec[(1+(nrow(TFBS_Z_Norm)*6)):(nrow(TFBS_Z_Norm)*7),2]
TFBS_Z_Norm[,9] = TFBS_DESeq2_Norm_vec[(1+(nrow(TFBS_Z_Norm)*7)):(nrow(TFBS_Z_Norm)*8),2]

ATAC_UP_TFBS_UP_TF_ATAC_HM = inner_join(ATAC_UP_TFBS_UP_TF,ATAC_Z_Norm,by="gene")
ATAC_UP_TFBS_DOWN_TF_ATAC_HM = inner_join(ATAC_UP_TFBS_DOWN_TF,ATAC_Z_Norm,by="gene")
ATAC_STEADY_TFBS_UP_TF_ATAC_HM = inner_join(ATAC_STEADY_TFBS_UP_TF,ATAC_Z_Norm,by="gene")
ATAC_STEADY_TFBS_DOWN_TF_ATAC_HM = inner_join(ATAC_STEADY_TFBS_DOWN_TF,ATAC_Z_Norm,by="gene")

ATAC_UP_TFBS_UP_TF_TFBS_HM = inner_join(ATAC_UP_TFBS_UP_TF,TFBS_Z_Norm,by="gene")
ATAC_UP_TFBS_DOWN_TF_TFBS_HM = inner_join(ATAC_UP_TFBS_DOWN_TF,TFBS_Z_Norm,by="gene")
ATAC_STEADY_TFBS_UP_TF_TFBS_HM = inner_join(ATAC_STEADY_TFBS_UP_TF,TFBS_Z_Norm,by="gene")
ATAC_STEADY_TFBS_DOWN_TF_TFBS_HM = inner_join(ATAC_STEADY_TFBS_DOWN_TF,TFBS_Z_Norm,by="gene")

HMdata_ATAC = rbind(ATAC_UP_TFBS_UP_TF_ATAC_HM,
               ATAC_UP_TFBS_DOWN_TF_ATAC_HM,
               ATAC_STEADY_TFBS_UP_TF_ATAC_HM,
               ATAC_STEADY_TFBS_DOWN_TF_ATAC_HM)

HMdata_TFBS = rbind(ATAC_UP_TFBS_UP_TF_TFBS_HM,
               ATAC_UP_TFBS_DOWN_TF_TFBS_HM,
               ATAC_STEADY_TFBS_UP_TF_TFBS_HM,
               ATAC_STEADY_TFBS_DOWN_TF_TFBS_HM)

options(repr.plot.width = 4, repr.plot.height = 10, repr.plot.res = 1000, repr.plot.pointsize = 10)
conte=HMdata_TFBS[,-1]
namen = HMdata_TFBS[,1]
rownames(conte)=as.vector(namen)
PCA=conte
pheatmap(PCA,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = F,
         legend = FALSE,
         border_color =NA,
         breaks=seq(from=-2,to=2,length.out = 100))

options(repr.plot.width = 4, repr.plot.height = 10, repr.plot.res = 1000, repr.plot.pointsize = 10)
conte=HMdata_ATAC[,-1]
namen = HMdata_ATAC[,1]
rownames(conte)=as.vector(namen)
PCA=conte
pheatmap(PCA,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = F,
         legend = FALSE,
         border_color =NA,
         breaks=seq(from=-2,to=2,length.out = 100))

