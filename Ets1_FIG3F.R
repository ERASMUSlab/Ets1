options(repr.plot.width = 3, repr.plot.height = 3, repr.plot.res = 1000, repr.plot.pointsize = 20)

Compartment_HistoneLevel = read.table("/data3/psg/Ets1/ProcessedData/Compartment_HistoneLevel.txt",sep = "\t", header=F,fill = TRUE)
violin_plot = Compartment_HistoneLevel[,c(13,7)]
colnames(violin_plot) = c("comp","value")

violin_plot = violin_plot[order(-violin_plot[,2]),]
violin_plot = violin_plot[(round(nrow(violin_plot)/10)):(round(nrow(violin_plot)/10)*9),]

ggplot(violin_plot, aes(x=factor(comp), y=value, fill=factor(comp)))+
geom_boxplot(fill=c("dodgerblue3","skyblue3","gold3","yellow3"),
             size=0.7,
             color='black',
             width=0.5,
             alpha=1)+
    guides(fill="none")+
theme_minimal()+
    labs(title="H3K27ac",
        subtitle = "")+
    xlab("")+
    ylab("Normalized active histone level\nat GM\n")+
theme(plot.title=element_text(face="bold.italic",hjust=0.5,size=15,color = "black"),
      axis.text.x=element_text(face="bold",size=12,color = "black"),
      axis.text.y=element_text(face="bold",size=10),
      plot.subtitle=element_text(face="bold",hjust=0.5,vjust=-3,size=0,color = "black"),
      axis.title.x = element_text(size=0,color = "black"),
      axis.title.y = element_text(face="bold",size=10,color = "black"))

options(repr.plot.width = 3, repr.plot.height = 3, repr.plot.res = 1000, repr.plot.pointsize = 20)

Compartment_HistoneLevel = read.table("/data3/psg/Ets1/ProcessedData/Compartment_HistoneLevel.txt",sep = "\t", header=F,fill = TRUE)
violin_plot = Compartment_HistoneLevel[,c(13,9)]
colnames(violin_plot) = c("comp","value")

violin_plot = violin_plot[order(-violin_plot[,2]),]
violin_plot = violin_plot[(round(nrow(violin_plot)/10)):(round(nrow(violin_plot)/10)*9),]

ggplot(violin_plot, aes(x=factor(comp), y=value, fill=factor(comp)))+
geom_boxplot(fill=c("dodgerblue3","skyblue3","gold3","yellow3"),
             size=0.7,
             color='black',
             width=0.5,
             alpha=1)+
    guides(fill="none")+
theme_minimal()+
    labs(title="H3K4me1",
        subtitle = "")+
    xlab("")+
    ylab("Normalized active histone level\nat GM\n")+
theme(plot.title=element_text(face="bold.italic",hjust=0.5,size=15,color = "black"),
      axis.text.x=element_text(face="bold",size=12,color = "black"),
      axis.text.y=element_text(face="bold",size=10),
      plot.subtitle=element_text(face="bold",hjust=0.5,vjust=-3,size=0,color = "black"),
      axis.title.x = element_text(size=0,color = "black"),
      axis.title.y = element_text(face="bold",size=10,color = "black"))

options(repr.plot.width = 3, repr.plot.height = 3, repr.plot.res = 1000, repr.plot.pointsize = 20)

Compartment_HistoneLevel = read.table("/data3/psg/Ets1/ProcessedData/Compartment_HistoneLevel.txt",sep = "\t", header=F,fill = TRUE)
violin_plot = Compartment_HistoneLevel[,c(13,11)]
colnames(violin_plot) = c("comp","value")

violin_plot = violin_plot[order(-violin_plot[,2]),]
violin_plot = violin_plot[(round(nrow(violin_plot)/10)):(round(nrow(violin_plot)/10)*9),]

ggplot(violin_plot, aes(x=factor(comp), y=value, fill=factor(comp)))+
geom_boxplot(fill=c("dodgerblue3","skyblue3","gold3","yellow3"),
             size=0.7,
             color='black',
             width=0.5,
             alpha=1)+
    guides(fill="none")+
theme_minimal()+
    labs(title="H3K4me3",
        subtitle = "")+
    xlab("")+
    ylab("Normalized active histone level\nat GM\n")+
theme(plot.title=element_text(face="bold.italic",hjust=0.5,size=15,color = "black"),
      axis.text.x=element_text(face="bold",size=12,color = "black"),
      axis.text.y=element_text(face="bold",size=10),
      plot.subtitle=element_text(face="bold",hjust=0.5,vjust=-3,size=0,color = "black"),
      axis.title.x = element_text(size=0,color = "black"),
      axis.title.y = element_text(face="bold",size=10,color = "black"))

