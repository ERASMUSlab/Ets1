BIN_DATA = read.table("/data3/psg/Ets1/ProcessedData/Ets1_p300_binDATA.txt",sep = "\t", header=T,fill = TRUE)

Qlevel = c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10")
Elevel = c("E1","E2","E3","E4","E5","E6","E7","E8","E9","E10")
label = c("A","B","C","D","E","F","G","H","I","J")

Qlevel = as.data.frame(Qlevel)
Elevel = as.data.frame(Elevel)
label = as.data.frame(label)


i = 1
Ets1_Q10 = BIN_DATA[which(BIN_DATA$Ets1_Q == Qlevel[i,1]),]
P1 = round(nrow(Ets1_Q10[which(Ets1_Q10$p300_Q == "Q1"),])/nrow(Ets1_Q10)*100)
P2 = round(nrow(Ets1_Q10[which(Ets1_Q10$p300_Q == "Q2"),])/nrow(Ets1_Q10)*100)
P3 = round(nrow(Ets1_Q10[which(Ets1_Q10$p300_Q == "Q3"),])/nrow(Ets1_Q10)*100)
P4 = round(nrow(Ets1_Q10[which(Ets1_Q10$p300_Q == "Q4"),])/nrow(Ets1_Q10)*100)
P5 = round(nrow(Ets1_Q10[which(Ets1_Q10$p300_Q == "Q5"),])/nrow(Ets1_Q10)*100)
P6 = round(nrow(Ets1_Q10[which(Ets1_Q10$p300_Q == "Q6"),])/nrow(Ets1_Q10)*100)
P7 = round(nrow(Ets1_Q10[which(Ets1_Q10$p300_Q == "Q7"),])/nrow(Ets1_Q10)*100)
P8 = round(nrow(Ets1_Q10[which(Ets1_Q10$p300_Q == "Q8"),])/nrow(Ets1_Q10)*100)
P9 = round(nrow(Ets1_Q10[which(Ets1_Q10$p300_Q == "Q9"),])/nrow(Ets1_Q10)*100)
P10 = round(nrow(Ets1_Q10[which(Ets1_Q10$p300_Q == "Q10"),])/nrow(Ets1_Q10)*100)

Ets1_to_p300 = data.frame(category=c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"),
                          count=c(P1,P2,P3,P4,P5,P6,P7,P8,P9,P10))

Ets1_to_p300$group = Elevel[i,1]
Ets1_to_p300$label = label[i,1]

for(i in 2:nrow(Qlevel)){
    
    Ets1_Q10 = BIN_DATA[which(BIN_DATA$Ets1_Q == Qlevel[i,1]),]
    P1 = round(nrow(Ets1_Q10[which(Ets1_Q10$p300_Q == "Q1"),])/nrow(Ets1_Q10)*100)
    P2 = round(nrow(Ets1_Q10[which(Ets1_Q10$p300_Q == "Q2"),])/nrow(Ets1_Q10)*100)
    P3 = round(nrow(Ets1_Q10[which(Ets1_Q10$p300_Q == "Q3"),])/nrow(Ets1_Q10)*100)
    P4 = round(nrow(Ets1_Q10[which(Ets1_Q10$p300_Q == "Q4"),])/nrow(Ets1_Q10)*100)
    P5 = round(nrow(Ets1_Q10[which(Ets1_Q10$p300_Q == "Q5"),])/nrow(Ets1_Q10)*100)
    P6 = round(nrow(Ets1_Q10[which(Ets1_Q10$p300_Q == "Q6"),])/nrow(Ets1_Q10)*100)
    P7 = round(nrow(Ets1_Q10[which(Ets1_Q10$p300_Q == "Q7"),])/nrow(Ets1_Q10)*100)
    P8 = round(nrow(Ets1_Q10[which(Ets1_Q10$p300_Q == "Q8"),])/nrow(Ets1_Q10)*100)
    P9 = round(nrow(Ets1_Q10[which(Ets1_Q10$p300_Q == "Q9"),])/nrow(Ets1_Q10)*100)
    P10 = round(nrow(Ets1_Q10[which(Ets1_Q10$p300_Q == "Q10"),])/nrow(Ets1_Q10)*100)
    
    Ets1_to_p300_frag = data.frame(category=c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"),
                                   count=c(P1,P2,P3,P4,P5,P6,P7,P8,P9,P10))
    
    Ets1_to_p300_frag$group = Elevel[i,1]
    Ets1_to_p300_frag$label = label[i,1]

    Ets1_to_p300 = rbind(Ets1_to_p300,Ets1_to_p300_frag)
    } 

options(repr.plot.width = 5, repr.plot.height = 3, repr.plot.res = 1000, repr.plot.pointsize = 40)
ggplot(Ets1_to_p300,aes(label, count)) +
geom_line(data=subset(Ets1_to_p300,category=="P1"), color="#E0F2F1", size=1, alpha = 1, group = "A") +
geom_line(data=subset(Ets1_to_p300,category=="P2"), color="#B2DFDB", size=1, alpha = 1, group = "B") +
geom_line(data=subset(Ets1_to_p300,category=="P3"), color="#80CBC4", size=1, alpha = 1, group = "C") +
geom_line(data=subset(Ets1_to_p300,category=="P4"), color="#4DB6AC", size=1, alpha = 1, group = "D") +
geom_line(data=subset(Ets1_to_p300,category=="P5"), color="#26A69A", size=1, alpha = 1, group = "E") +
geom_line(data=subset(Ets1_to_p300,category=="P6"), color="#009688", size=1, alpha = 1, group = "F") +
geom_line(data=subset(Ets1_to_p300,category=="P7"), color="#00897B", size=1, alpha = 1, group = "G") +
geom_line(data=subset(Ets1_to_p300,category=="P8"), color="#00796B", size=1, alpha = 1, group = "H") +
geom_line(data=subset(Ets1_to_p300,category=="P9"), color="#00695C", size=1, alpha = 1, group = "I") +
geom_line(data=subset(Ets1_to_p300,category=="P10"), color="#004D40", size=1, alpha = 1, group = "J") +

geom_point(data=subset(Ets1_to_p300,category=="P1"), color="#E0F2F1", size=2, alpha = 1, group = "A",shape=21,fill="white") +
geom_point(data=subset(Ets1_to_p300,category=="P2"), color="#B2DFDB", size=2, alpha = 1, group = "B",shape=21,fill="white") +
geom_point(data=subset(Ets1_to_p300,category=="P3"), color="#80CBC4", size=2, alpha = 1, group = "C",shape=21,fill="white") +
geom_point(data=subset(Ets1_to_p300,category=="P4"), color="#4DB6AC", size=2, alpha = 1, group = "D",shape=21,fill="white") +
geom_point(data=subset(Ets1_to_p300,category=="P5"), color="#26A69A", size=2, alpha = 1, group = "E",shape=21,fill="white") +
geom_point(data=subset(Ets1_to_p300,category=="P6"), color="#009688", size=2, alpha = 1, group = "F",shape=21,fill="white") +
geom_point(data=subset(Ets1_to_p300,category=="P7"), color="#00897B", size=2, alpha = 1, group = "G",shape=21,fill="white") +
geom_point(data=subset(Ets1_to_p300,category=="P8"), color="#00796B", size=2, alpha = 1, group = "H",shape=21,fill="white") +
geom_point(data=subset(Ets1_to_p300,category=="P9"), color="#00695C", size=2, alpha = 1, group = "I",shape=21,fill="white") +
geom_point(data=subset(Ets1_to_p300,category=="P10"), color="#004D40", size=2, alpha = 1, group = "J",shape=21,fill="white") +
scale_x_discrete(labels = c("E1","E2","E3","E4","E5","E6","E7","E8","E9","E10"))+
theme_classic(base_size = 10) +
ggtitle(NULL)+
ylab("p300 quantile distribution\n(Ratio,%)") +
xlab("Ets1 binding level") +
theme(plot.title=element_text(face="bold",hjust=0.5,size=25,color = "black"),
      axis.text.x=element_text(face="bold",size=10),
      axis.text.y=element_text(face="bold",size=10),
      axis.title.x = element_text(face="bold",size = 13,color = "black"),
      axis.title.y = element_text(face="bold",size = 13,color = "black"),
      legend.title=element_text(face="bold",size=15), 
      legend.text=element_text(face="bold",size=15))

