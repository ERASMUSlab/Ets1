ATAC_Full_data = read.table("/Ets1/ProcessedData/Diff_Access_TF_Full_data.txt",sep = "\t", header=T,fill = TRUE)
TFBS_Full_data = read.table("/Ets1/ProcessedData/Diff_TFBS_TF_Full_data.txt",sep = "\t", header=T,fill = TRUE)
ATAC_DESeq2_Norm = read.table("/Ets1/ProcessedData/Access_TF_count_data.txt",sep = "\t", header=T,fill = TRUE)
TFBS_DESeq2_Norm = read.table("/Ets1/ProcessedData/TFBS_TF_count_data.txt",sep = "\t", header=T,fill = TRUE)
RAWDATA = read.table("/Ets1/ProcessedData/Stage1_RNAcount.txt",sep = "\t", header=T,fill = TRUE)

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

ATAC_UP = ATAC_Full_data[which(ATAC_Full_data$log2FoldChange > -log2(1.2)),]
ATAC_UP = ATAC_UP[which(ATAC_UP$pvalue < 0.05),]
ATAC_UP = ATAC_UP[order(-ATAC_UP$log2FoldChange),]


ATAC_UP_fit = ATAC_UP[,c(1,3)]
ATAC_UP_fit = as.data.frame(ATAC_UP_fit)
colnames(ATAC_UP_fit) = c("TF","ATAC_LFC")

TFBS_UP_fit = TFBS_UP[,c(1,3)]
TFBS_UP_fit = as.data.frame(TFBS_UP_fit)
colnames(TFBS_UP_fit) = c("TF","TFBS_LFC")

ATAC_UP_TFBS_UP_fit = inner_join(ATAC_UP_fit,TFBS_UP_fit,by="TF")

ATAC_UP_TFBS_UP_fit$ATAC_FC = 2^(ATAC_UP_TFBS_UP_fit[,2])
ATAC_UP_TFBS_UP_fit$TFBS_FC = 2^(ATAC_UP_TFBS_UP_fit[,3])

ATAC_UP_TFBS_UP_fit = ATAC_UP_TFBS_UP_fit[order(-ATAC_UP_TFBS_UP_fit$ATAC_LFC),]

RAWDATA_fit = RAWDATA[,1:9]

RAWDATA_fit$sum = RAWDATA_fit[,2]+RAWDATA_fit[,3]+RAWDATA_fit[,4]+RAWDATA_fit[,5]+RAWDATA_fit[,6]+RAWDATA_fit[,7]
RAWDATA_fit = RAWDATA_fit[which(RAWDATA_fit$sum > 10),]
RAWDATA_fit = RAWDATA_fit[,1]
RAWDATA_fit = as.data.frame(RAWDATA_fit)
colnames(RAWDATA_fit) = "TF"
ATAC_UP_TFBS_UP_fit = inner_join(ATAC_UP_TFBS_UP_fit,RAWDATA_fit,by="TF")

ATAC = as.data.frame(cbind(ATAC_UP_TFBS_UP_fit[1,1],"ATAC",ATAC_UP_TFBS_UP_fit[1,4]))
TFBS = as.data.frame(cbind(ATAC_UP_TFBS_UP_fit[1,1],"TFBS",ATAC_UP_TFBS_UP_fit[1,5]))

bar_DF = rbind(ATAC,TFBS)

for(i in 2:nrow(ATAC_UP_TFBS_UP_fit)){
ATAC = as.data.frame(cbind(ATAC_UP_TFBS_UP_fit[i,1],"ATAC",ATAC_UP_TFBS_UP_fit[i,4]))
TFBS = as.data.frame(cbind(ATAC_UP_TFBS_UP_fit[i,1],"TFBS",ATAC_UP_TFBS_UP_fit[i,5]))

bar_DF_frag = rbind(ATAC,TFBS)

    bar_DF = rbind(bar_DF,bar_DF_frag)
    
    }
colnames(bar_DF)=c("TF","cate","value")
bar_DF[,3] = as.numeric(bar_DF[,3])
bar_DF$RANK = c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11)
bar_DF[,4] = as.factor(bar_DF[,4])
label = unique(as.data.frame(bar_DF[,1]))

options(repr.plot.width = 8, repr.plot.height = 3, repr.plot.res = 1000, repr.plot.pointsize = 40)

ggplot(bar_DF,aes(x=RANK,y=value,fill=cate))+
geom_bar(stat="identity",position="dodge",width=0.6)+
scale_fill_manual(values = c("darkgreen","purple1"))+
theme_classic(base_size = 10) +
scale_x_discrete(labels = label[,1])+
ggtitle(label = "",
        subtitle = "") +
ylab("Chromatin accessibility FC\n") +
theme(plot.title=element_text(face="bold.italic",hjust=0.5,size=20,color = "black"),
      plot.subtitle=element_text(face="bold.italic",hjust=0.5,size=10,color = "black"),
      axis.text.x=element_text(face="bold",size=10,color = "black"),
      axis.text.y=element_text(face="bold",size=12,color = "black"),
      axis.title.x = element_text(face="bold",size = 0,color = "darkblue"),
      axis.title.y = element_text(face="bold",size = 12,color = "black"),
      legend.title=element_text(face="bold",size=10), 
      legend.text=element_text(size=8),
      legend.position = "none")

