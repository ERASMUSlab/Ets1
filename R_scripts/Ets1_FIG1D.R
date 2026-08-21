Diff_0TO4_DEG_UP_GO = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG1/3_DEGs/DEG_result/Diff_0TO4_DEG_UP_GO.txt",sep = "\t", header=T,fill = TRUE)
Diff_4TO7_DEG_UP_GO = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG1/3_DEGs/DEG_result/Diff_4TO7_DEG_UP_GO.txt",sep = "\t", header=T,fill = TRUE)
Diff_7TO10_DEG_UP_GO = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG1/3_DEGs/DEG_result/Diff_7TO10_DEG_UP_GO.txt",sep = "\t", header=T,fill = TRUE)

Diff_0TO4_DEG_UP_GO = Diff_0TO4_DEG_UP_GO[,c(2,10)]
Diff_4TO7_DEG_UP_GO = Diff_4TO7_DEG_UP_GO[,c(2,10)]
Diff_7TO10_DEG_UP_GO = Diff_7TO10_DEG_UP_GO[,c(2,10)]

colnames(Diff_0TO4_DEG_UP_GO)[2] = "RF_0TO4"
colnames(Diff_4TO7_DEG_UP_GO)[2] = "RF_4TO7"
colnames(Diff_7TO10_DEG_UP_GO)[2] = "RF_7TO10"

Diff_DEG_UP_GO = full_join(Diff_0TO4_DEG_UP_GO,Diff_4TO7_DEG_UP_GO,by="Description")
Diff_DEG_UP_GO = full_join(Diff_DEG_UP_GO,Diff_7TO10_DEG_UP_GO,by="Description")

Diff_DEG_UP_GO[is.na(Diff_DEG_UP_GO)] = 0

Diff_DEG_UP_GO$label = "NOT"

for(i in 1:nrow(Diff_DEG_UP_GO)){
    RF_0TO4 = Diff_DEG_UP_GO[i,2]
    RF_4TO7 = Diff_DEG_UP_GO[i,3]
    RF_7TO10 = Diff_DEG_UP_GO[i,4]

    if(RF_0TO4 >= RF_4TO7){
        if(RF_0TO4 >= RF_7TO10){
             Diff_DEG_UP_GO[i,5] = "RF_0TO4"
            }
        }

    if(RF_4TO7 >= RF_0TO4){
        if(RF_4TO7 >= RF_7TO10){
             Diff_DEG_UP_GO[i,5] = "RF_4TO7"
            }
        }

    if(RF_7TO10 >= RF_0TO4){
        if(RF_7TO10 >= RF_4TO7){
             Diff_DEG_UP_GO[i,5] = "RF_7TO10"
            }
        }
    }


Diff_DEG_UP_GO_RF_0TO4 = Diff_DEG_UP_GO[which(Diff_DEG_UP_GO[,5] == "RF_0TO4"),]
Diff_DEG_UP_GO_RF_4TO7 = Diff_DEG_UP_GO[which(Diff_DEG_UP_GO[,5] == "RF_4TO7"),]
Diff_DEG_UP_GO_RF_7TO10 = Diff_DEG_UP_GO[which(Diff_DEG_UP_GO[,5] == "RF_7TO10"),]

Diff_DEG_UP_GO_plot = rbind(Diff_DEG_UP_GO_RF_0TO4[which(Diff_DEG_UP_GO_RF_0TO4[,1] == "ossification"),],
                            Diff_DEG_UP_GO_RF_0TO4[which(Diff_DEG_UP_GO_RF_0TO4[,1] == "bone development"),],
                            Diff_DEG_UP_GO_RF_0TO4[which(Diff_DEG_UP_GO_RF_0TO4[,1] == "osteoblast differentiation"),],
                            Diff_DEG_UP_GO_RF_4TO7[which(Diff_DEG_UP_GO_RF_4TO7[,1] == "canonical glycolysis"),],
                            Diff_DEG_UP_GO_RF_4TO7[which(Diff_DEG_UP_GO_RF_4TO7[,1] == "response to vitamin"),],
                            Diff_DEG_UP_GO_RF_4TO7[which(Diff_DEG_UP_GO_RF_4TO7[,1] == "glucose catabolic process"),],
                            Diff_DEG_UP_GO_RF_7TO10[which(Diff_DEG_UP_GO_RF_7TO10[,1] == "cell adhesion mediated by integrin"),],
                            Diff_DEG_UP_GO_RF_7TO10[which(Diff_DEG_UP_GO_RF_7TO10[,1] == "regulation of cell-substrate adhesion"),],
                            Diff_DEG_UP_GO_RF_7TO10[which(Diff_DEG_UP_GO_RF_7TO10[,1] == "regulation of collagen biosynthetic process"),])

plot_DF = rbind(cbind(paste(9,Diff_DEG_UP_GO_plot[1,1],sep= "  "),"3",Diff_DEG_UP_GO_plot[1,2]),
                cbind(paste(9,Diff_DEG_UP_GO_plot[1,1],sep= "  "),"2",Diff_DEG_UP_GO_plot[1,3]),
                cbind(paste(9,Diff_DEG_UP_GO_plot[1,1],sep= "  "),"1",Diff_DEG_UP_GO_plot[1,4]))

for(i in 2:nrow(Diff_DEG_UP_GO_plot)){
    plot_DF_frag = rbind(cbind(paste(10-i,Diff_DEG_UP_GO_plot[i,1],sep= "  "),"3",Diff_DEG_UP_GO_plot[i,2]),
                         cbind(paste(10-i,Diff_DEG_UP_GO_plot[i,1],sep= "  "),"2",Diff_DEG_UP_GO_plot[i,3]),
                         cbind(paste(10-i,Diff_DEG_UP_GO_plot[i,1],sep= "  "),"1",Diff_DEG_UP_GO_plot[i,4]))

    plot_DF = rbind(plot_DF,plot_DF_frag)
    }

plot_DF = as.data.frame(plot_DF)
colnames(plot_DF) = c("Description","cate","value")
plot_DF[,3] = as.numeric(plot_DF[,3])

options(repr.plot.width = 12, repr.plot.height = 6, repr.plot.res = 1000, repr.plot.pointsize =1)

ggplot(plot_DF,aes(value, Description, fill=cate)) +
geom_bar(stat="identity",position="dodge",colour="black",width=0.85,linewidth=0.1)+
theme_dose(0) +
ylab("") +
xlab("Enrichment score") +
scale_y_discrete(labels = c("regulation of collagen biosynthetic process", 
                            "regulation of cell-substrate adhesion", 
                            "cell adhesion mediated by integrin",
                            "glucose catabolic process", 
                            "response to vitamin", 
                            "canonical glycolysis",
                            "osteoblast differentiation", 
                            "bone development", 
                            "ossification"))+
scale_fill_manual(values = c("skyblue3","tan3","purple"))+
theme(plot.title=element_text(face="bold.italic",hjust=0.5,size=20,color = "black"),
      plot.subtitle=element_text(face="bold.italic",hjust=0.5,size=10,color = "black"),
      axis.text.x=element_text(face="bold",size=15,color = "black"),
      axis.text.y=element_text(face="bold",size=20,color = "black"),
      axis.title.x = element_text(face="bold",size = 20,color = "black"),
      axis.title.y = element_text(face="bold",size = 20,color = "black"),
      legend.title=element_text(face="bold",size=10), 
      legend.text=element_text(face="bold",size = 9,color = "black"),
      legend.position = "none") 

