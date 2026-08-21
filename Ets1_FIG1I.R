mm10_RE_ChIPRNA = read.table("/data3/psg/NGS_2025/4DN/HiC_data/FiG1/DBs/mm10_RE_ChIPRNA.txt", sep = "\t", header=T,fill = TRUE)

mm10_RE_ChIPRNA$H3K27ac_set = "red"

for(i in 1:nrow(mm10_RE_ChIPRNA)){
    
    if(mm10_RE_ChIPRNA[i,11] >= 1.5){
        if(mm10_RE_ChIPRNA[i,14] >= 1.5){
            mm10_RE_ChIPRNA[i,16] = "green"
            }
        }

    if(mm10_RE_ChIPRNA[i,11] >= 1.5){
        if(mm10_RE_ChIPRNA[i,14] <= -1.5){
            mm10_RE_ChIPRNA[i,16] = "red"
            }
        }

    if(mm10_RE_ChIPRNA[i,11] <= -1.5){
        if(mm10_RE_ChIPRNA[i,14] >= 1.5){
            mm10_RE_ChIPRNA[i,16] = "red"
            }
        }

    if(mm10_RE_ChIPRNA[i,11] <= -1.5){
        if(mm10_RE_ChIPRNA[i,14] <= -1.5){
            mm10_RE_ChIPRNA[i,16] = "green"
            }
        }

    if(mm10_RE_ChIPRNA[i,11] >= -1.5){
        if(mm10_RE_ChIPRNA[i,11] <= 1.5){
            if(mm10_RE_ChIPRNA[i,14] >= -1.5){
                if(mm10_RE_ChIPRNA[i,14] <= 1.5){
                    mm10_RE_ChIPRNA[i,16] = "green"
                    }
                }
            }
        }
    }

mm10_RE_ChIPRNA$H3K4me1_set = "red"

for(i in 1:nrow(mm10_RE_ChIPRNA)){
    
    if(mm10_RE_ChIPRNA[i,12] >= 1.5){
        if(mm10_RE_ChIPRNA[i,14] >= 1.5){
            mm10_RE_ChIPRNA[i,17] = "green"
            }
        }

    if(mm10_RE_ChIPRNA[i,12] >= 1.5){
        if(mm10_RE_ChIPRNA[i,14] <= -1.5){
            mm10_RE_ChIPRNA[i,17] = "red"
            }
        }

    if(mm10_RE_ChIPRNA[i,12] <= -1.5){
        if(mm10_RE_ChIPRNA[i,14] >= 1.5){
            mm10_RE_ChIPRNA[i,17] = "red"
            }
        }

    if(mm10_RE_ChIPRNA[i,12] <= -1.5){
        if(mm10_RE_ChIPRNA[i,14] <= -1.5){
            mm10_RE_ChIPRNA[i,17] = "green"
            }
        }

    if(mm10_RE_ChIPRNA[i,12] >= -1.5){
        if(mm10_RE_ChIPRNA[i,12] <= 1.5){
            if(mm10_RE_ChIPRNA[i,14] >= -1.5){
                if(mm10_RE_ChIPRNA[i,14] <= 1.5){
                    mm10_RE_ChIPRNA[i,17] = "green"
                    }
                }
            }
        }
    }

mm10_RE_ChIPRNA$H3K4me3_set = "red"

for(i in 1:nrow(mm10_RE_ChIPRNA)){
    
    if(mm10_RE_ChIPRNA[i,13] >= 1.5){
        if(mm10_RE_ChIPRNA[i,14] >= 1.5){
            mm10_RE_ChIPRNA[i,18] = "green"
            }
        }

    if(mm10_RE_ChIPRNA[i,13] >= 1.5){
        if(mm10_RE_ChIPRNA[i,14] <= -1.5){
            mm10_RE_ChIPRNA[i,18] = "red"
            }
        }

    if(mm10_RE_ChIPRNA[i,13] <= -1.5){
        if(mm10_RE_ChIPRNA[i,14] >= 1.5){
            mm10_RE_ChIPRNA[i,18] = "red"
            }
        }

    if(mm10_RE_ChIPRNA[i,13] <= -1.5){
        if(mm10_RE_ChIPRNA[i,14] <= -1.5){
            mm10_RE_ChIPRNA[i,18] = "green"
            }
        }

    if(mm10_RE_ChIPRNA[i,13] >= -1.5){
        if(mm10_RE_ChIPRNA[i,13] <= 1.5){
            if(mm10_RE_ChIPRNA[i,14] >= -1.5){
                if(mm10_RE_ChIPRNA[i,14] <= 1.5){
                    mm10_RE_ChIPRNA[i,18] = "green"
                    }
                }
            }
        }
    }

for(i in 1:nrow(mm10_RE_ChIPRNA)){
    if(mm10_RE_ChIPRNA[i,11]>0){mm10_RE_ChIPRNA[i,11] = mm10_RE_ChIPRNA[i,11]-1}
    if(mm10_RE_ChIPRNA[i,11]<0){mm10_RE_ChIPRNA[i,11] = mm10_RE_ChIPRNA[i,11]+1}

    if(mm10_RE_ChIPRNA[i,12]>0){mm10_RE_ChIPRNA[i,12] = mm10_RE_ChIPRNA[i,12]-1}
    if(mm10_RE_ChIPRNA[i,12]<0){mm10_RE_ChIPRNA[i,12] = mm10_RE_ChIPRNA[i,12]+1}

    if(mm10_RE_ChIPRNA[i,13]>0){mm10_RE_ChIPRNA[i,13] = mm10_RE_ChIPRNA[i,13]-1}
    if(mm10_RE_ChIPRNA[i,13]<0){mm10_RE_ChIPRNA[i,13] = mm10_RE_ChIPRNA[i,13]+1}

    if(mm10_RE_ChIPRNA[i,14]>0){mm10_RE_ChIPRNA[i,14] = mm10_RE_ChIPRNA[i,14]-1}
    if(mm10_RE_ChIPRNA[i,14]<0){mm10_RE_ChIPRNA[i,14] = mm10_RE_ChIPRNA[i,14]+1}
    }

round(nrow(mm10_RE_ChIPRNA[which(mm10_RE_ChIPRNA$H3K27ac_set == "green"),])/nrow(mm10_RE_ChIPRNA)*100,2)
round(nrow(mm10_RE_ChIPRNA[which(mm10_RE_ChIPRNA$H3K27ac_set == "red"),])/nrow(mm10_RE_ChIPRNA)*100,2)

round(nrow(mm10_RE_ChIPRNA[which(mm10_RE_ChIPRNA$H3K4me1_set == "green"),])/nrow(mm10_RE_ChIPRNA)*100,2)
round(nrow(mm10_RE_ChIPRNA[which(mm10_RE_ChIPRNA$H3K4me1_set == "red"),])/nrow(mm10_RE_ChIPRNA)*100,2)

round(nrow(mm10_RE_ChIPRNA[which(mm10_RE_ChIPRNA$H3K4me3_set == "green"),])/nrow(mm10_RE_ChIPRNA)*100,2)
round(nrow(mm10_RE_ChIPRNA[which(mm10_RE_ChIPRNA$H3K4me3_set == "red"),])/nrow(mm10_RE_ChIPRNA)*100,2)

options(repr.plot.width = 8, repr.plot.height = 8, repr.plot.res = 1000, repr.plot.pointsize = 10)

H3K27acFIG = ggplot(mm10_RE_ChIPRNA, aes(x=H3K27ac_FC,y=RNA_FC)) + 
geom_point(data=subset(mm10_RE_ChIPRNA,H3K27ac_set=="red"), color="darkred", size=1) +
geom_point(data=subset(mm10_RE_ChIPRNA,H3K27ac_set=="green"), color="darkgreen", size=1) +
geom_hline(yintercept=c(-0.5,0.5),
                         color = "black",
                         size = 0.5) +
geom_vline(xintercept=c(-0.5,0.5),
                         color = "black",
                         size = 0.5) +
ylim(-10,10) +
xlim(-10,10) +
theme_classic(base_size = 10) + 
theme(legend.position = "none")+
xlab("H3K27ac FC") +
ylab("\nExpression FC") +
theme(plot.title=element_text(face="bold",hjust=0.5,size=15,color = "black"),
      axis.text.x=element_text(face="bold",size=15,color = "black"),
      axis.text.y=element_text(face="bold",size=15,color = "black"),
      axis.title.x = element_text(face="bold",size = 20,color = "black"),
      axis.title.y = element_text(face="bold",size = 20,color = "black"),
      legend.title=element_text(face="bold",size=10), 
      legend.text=element_text(face="bold",size=8))

H3K4me1FIG = ggplot(mm10_RE_ChIPRNA, aes(x=H3K4me1_FC,y=RNA_FC)) + 
geom_point(data=subset(mm10_RE_ChIPRNA,H3K4me1_set=="red"), color="darkred", size=1) +
geom_point(data=subset(mm10_RE_ChIPRNA,H3K4me1_set=="green"), color="darkgreen", size=1) +
geom_hline(yintercept=c(-0.5,0.5),
                         color = "black",
                         size = 0.5) +
geom_vline(xintercept=c(-0.5,0.5),
                         color = "black",
                         size = 0.5) +
ylim(-10,10) +
xlim(-10,10) +
theme_classic(base_size = 10) + 
theme(legend.position = "none")+
xlab("H3K4me1 FC") +
ylab("\nExpression FC") +
theme(plot.title=element_text(face="bold",hjust=0.5,size=15,color = "black"),
      axis.text.x=element_text(face="bold",size=15,color = "black"),
      axis.text.y=element_text(face="bold",size=15,color = "black"),
      axis.title.x = element_text(face="bold",size = 20,color = "black"),
      axis.title.y = element_text(face="bold",size = 20,color = "black"),
      legend.title=element_text(face="bold",size=10), 
      legend.text=element_text(face="bold",size=8))

H3K4me3FIG = ggplot(mm10_RE_ChIPRNA, aes(x=H3K4me3_FC,y=RNA_FC)) + 
geom_point(data=subset(mm10_RE_ChIPRNA,H3K4me3_set=="red"), color="darkred", size=1) +
geom_point(data=subset(mm10_RE_ChIPRNA,H3K4me3_set=="green"), color="darkgreen", size=1) +
geom_hline(yintercept=c(-0.5,0.5),
                         color = "black",
                         size = 0.5) +
geom_vline(xintercept=c(-0.5,0.5),
                         color = "black",
                         size = 0.5) +
ylim(-10,10) +
xlim(-10,10) +
theme_classic(base_size = 10) + 
theme(legend.position = "none")+
xlab("H3K4me3 FC") +
ylab("\nExpression FC") +
theme(plot.title=element_text(face="bold",hjust=0.5,size=15,color = "black"),
      axis.text.x=element_text(face="bold",size=15,color = "black"),
      axis.text.y=element_text(face="bold",size=15,color = "black"),
      axis.title.x = element_text(face="bold",size = 20,color = "black"),
      axis.title.y = element_text(face="bold",size = 20,color = "black"),
      legend.title=element_text(face="bold",size=10), 
      legend.text=element_text(face="bold",size=8))

ppp=arrangeGrob(H3K27acFIG,
                H3K4me1FIG,
                H3K4me3FIG,

    ncol = 3,
    nrow = 1)

options(repr.plot.width = 16, repr.plot.height = 5, repr.plot.res = 300, repr.plot.pointsize = 40)
grid.arrange(ppp)

