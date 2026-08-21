RAWDATA = read.table("/Ets1/ProcessedData/Bulk_RNAseq_expression_counts.txt",sep = "\t", header=T,fill = TRUE)

condition = c("DAY0","DAY0","DAY0","DAY0",
              "DAY4","DAY4","DAY4","DAY4",
              "DAY7","DAY7","DAY7","DAY7",
              "DAY10","DAY10","DAY10","DAY10")

label = c("DAY0","","","",
          "","DAY4","","",
          "","","DAY7","",
          "","","DAY10","")

predata = t(RAWDATA)
predata=as.data.frame(predata)
predata_contents = predata
predata_contents = predata_contents
predata_contents = predata_contents[-1,]
predata_colname = as.vector(predata[1,])
colnames(predata_contents)=predata_colname
predata=predata_contents
data=as.data.frame(predata)

for(i in 1:ncol(data)){
    data[,i]=as.numeric(data[,i])
}

datapre = data
dataset=datapre[ , which(apply(datapre, 2, var) != 0)]
pca <- prcomp(dataset, center = TRUE,scale. = TRUE)
pca_perc=round(100*pca$sdev^2/sum(pca$sdev^2),1)


df_pca_data=data.frame(PC1 = pca$x[,1], 
                       PC2 = pca$x[,2], 
                       PC3 = pca$x[,3], 
                       sample = rownames(dataset), 
                       condition=condition)

options(repr.plot.width = 5, repr.plot.height = 5, repr.plot.res = 1000, repr.plot.pointsize = 20)

ggplot(df_pca_data,aes(PC1,PC2,color = condition))+
geom_point(size=5)+
labs(x=paste0("PC1 (",pca_perc[1],")"), 
     y=paste0("PC2 (",pca_perc[2],")"))+
ylim(-130,130)+
xlim(-130,130)+
theme_classic(base_size = 10) +
scale_color_manual(values = c("pink3","purple","skyblue3","tan3")) +
theme(plot.title=element_text(face="bold",hjust=0.5,size=25,color = "darkblue"),
      axis.text.x=element_text(face="bold",size=12),
      axis.text.y=element_text(face="bold",size=12),
      axis.title.x = element_text(face="bold",size = 15),
      axis.title.y = element_text(face="bold",size = 15),
      legend.title=element_text(face="bold",size=0,color = "darkblue"),
      legend.text=element_text(face="bold",size=8))+
theme_classic(base_size = 10) +
geom_text(aes(label=label),hjust=0.5, vjust=-1,size=6, fontface = "bold")+
guides(color="none")

msigdbr = msigdbr(species = "mouse", category = "C5", subcategory = "GO:BP")[,3:4]
OSSI_BIG_SET = msigdbr[grepl("OSSI",msigdbr$gs_name), ]

OSSI_BIG_SET_term = OSSI_BIG_SET[,1]
OSSI_BIG_SET_term = as.data.frame(OSSI_BIG_SET_term)
colnames(OSSI_BIG_SET_term) = "Term"
OSSI_BIG_SET_term = unique(OSSI_BIG_SET_term)
OSSI_BIG_SET_term$GENE_NUM = 0

for(i in 1:nrow(OSSI_BIG_SET_term)){
    OSSI_BIG_SET_term[i,2] = nrow(unique(OSSI_BIG_SET[which(OSSI_BIG_SET[,1] == OSSI_BIG_SET_term[i,1]),]))
    }

OSSI_GOBP_term = OSSI_BIG_SET_term[order(-OSSI_BIG_SET_term$GENE_NUM),]
OSSI_GOBP_DB = OSSI_BIG_SET
OSSI_GOBP_gene = OSSI_GOBP_DB[,2]
OSSI_GOBP_gene = as.data.frame(OSSI_GOBP_gene)
OSSI_GOBP_gene = unique(OSSI_GOBP_gene)
colnames(OSSI_GOBP_gene) = "gene"

RAWDATA = read.table("/Ets1/ProcessedData/Bulk_RNAseq_expression_counts.txt",sep = "\t", header=T,fill = TRUE)

RAWDATA = inner_join(RAWDATA,OSSI_GOBP_gene,by="gene")

condition = c("DAY0","DAY0","DAY0","DAY0",
              "DAY4","DAY4","DAY4","DAY4",
              "DAY7","DAY7","DAY7","DAY7",
              "DAY10","DAY10","DAY10","DAY10")

label = c("DAY0","","","",
          "","DAY4","","",
          "","","DAY7","",
          "","","","DAY10")
    
predata = t(RAWDATA)
predata=as.data.frame(predata)
predata_contents = predata
predata_contents = predata_contents
predata_contents = predata_contents[-1,]
predata_colname = as.vector(predata[1,])
colnames(predata_contents)=predata_colname
predata=predata_contents
data=as.data.frame(predata)

for(i in 1:ncol(data)){
    data[,i]=as.numeric(data[,i])
}

datapre = data
dataset=datapre[ , which(apply(datapre, 2, var) != 0)]
pca <- prcomp(dataset, center = TRUE,scale. = TRUE)
pca_perc=round(100*pca$sdev^2/sum(pca$sdev^2),1)


df_pca_data=data.frame(PC1 = pca$x[,1], 
                       PC2 = pca$x[,2], 
                       PC3 = pca$x[,3], 
                       sample = rownames(dataset), 
                       condition=condition)

options(repr.plot.width = 5, repr.plot.height = 5, repr.plot.res = 1000, repr.plot.pointsize = 20)

ggplot(df_pca_data,aes(PC1,PC2,color = condition))+
geom_point(size=5)+
labs(x=paste0("PC1 (",pca_perc[1],")"), 
     y=paste0("PC2 (",pca_perc[2],")"))+
ylim(-30,30)+
xlim(-30,30)+
theme_classic(base_size = 10) +
scale_color_manual(values = c("pink3","purple","skyblue3","tan3")) +
theme(plot.title=element_text(face="bold",hjust=0.5,size=25,color = "darkblue"),
      axis.text.x=element_text(face="bold",size=12),
      axis.text.y=element_text(face="bold",size=12),
      axis.title.x = element_text(face="bold",size = 15),
      axis.title.y = element_text(face="bold",size = 15),
      legend.title=element_text(face="bold",size=0,color = "darkblue"),
      legend.text=element_text(face="bold",size=8))+
theme_classic(base_size = 10) +
geom_text(aes(label=label),hjust=0.5, vjust=-1,size=6, fontface = "bold")+
guides(color="none")

