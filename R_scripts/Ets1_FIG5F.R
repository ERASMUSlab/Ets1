UMAP_DIFF_correlation_50_2_1_1 = as.data.frame(vroom("/Ets1/ProcessedData/UMAP_DIFF_correlation_50_2_1_1.txt", delim = "\t",show_col_types = FALSE))
UMAP_DIFF_correlation_50_2_1_1_fit = UMAP_DIFF_correlation_50_2_1_1[,1:2]
UMAP_DIFF_correlation_50_2_1_1_label = paste(UMAP_DIFF_correlation_50_2_1_1[,3],UMAP_DIFF_correlation_50_2_1_1[,4],UMAP_DIFF_correlation_50_2_1_1[,5])
UMAP_DIFF_correlation_50_2_1_1_label = as.data.frame(UMAP_DIFF_correlation_50_2_1_1_label)
colnames(UMAP_DIFF_correlation_50_2_1_1_label) = "label"

result <- kmeans(UMAP_DIFF_correlation_50_2_1_1_fit, centers = 22, iter.max = 10000, algorithm = "MacQueen")

UMAP_DIFF_correlation_50_2_1_1 = cbind(UMAP_DIFF_correlation_50_2_1_1,as.data.frame(result$cluster))
colnames(UMAP_DIFF_correlation_50_2_1_1)[ncol(UMAP_DIFF_correlation_50_2_1_1)] = "cluster"

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

colnames(Q1_data) = c("cluster","Q_ratio","label","sort")
colnames(Q10_data) = c("cluster","Q_ratio","label","sort")
QQ_data = rbind(Q1_data,Q10_data)

options(repr.plot.width = 5, repr.plot.height = 3, repr.plot.res = 1000, repr.plot.pointsize = 40)
ggplot(QQ_data,aes(sort, y = Q_ratio,)) +
geom_line(data=subset(QQ_data,label=="Q1"), color="darkblue", size=0.7, alpha = 1, group = "Q1") +
geom_line(data=subset(QQ_data,label=="Q10"), color="darkred", size=0.7, alpha = 1, group = "Q10") +
geom_point(color = "black", size = 2) + 
theme_classic(base_size = 10) +
scale_x_discrete(labels = QQ_data[,1])+ 
ggtitle(NULL)+
xlab("Cluster") +
ylab("Ratio of bins\n(Q10,Q1)") +
theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 25, color = "black"),
      axis.text.x = element_text(face = "bold", size = 8),
      axis.text.y = element_text(face = "bold", size = 8),
      axis.title.x = element_text(face = "bold", size = 10, color = "black"),
      axis.title.y = element_text(face = "bold", size = 10, color = "black"),
      legend.position = "none",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(face = "bold", size = 7))

