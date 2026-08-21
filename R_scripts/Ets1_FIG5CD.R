UMAP_function= function(metric_para,n_neighbors_para,n_components_para,min_dist_para,spread_para,n_threads_para){

MAT_forUMAP = read.table("/data3/psg/Ets1/ProcessedData/MATRIX_DIFF.txt",sep = "\t", header=T)
MATRIX_DAY0_NOR_filltered2 = read.table("/data3/psg/Ets1/ProcessedData/BIN_AB_Ets1_p300_data.txt",sep = "\t", header=F)

UMAP <- uwot::umap(MAT_forUMAP,
                   metric = metric_para,
                   n_neighbors = n_neighbors_para,
                   n_components = n_components_para,
                   min_dist = min_dist_para,
                   spread = spread_para,
                   n_threads = n_threads_para,
                   seed = 7777,
                   verbose = FALSE)

UMAP = as.data.frame(UMAP)
UMAP = UMAP[,1:2]
colnames(UMAP) = c("X","Y")
UMAP = cbind(UMAP,MATRIX_DAY0_NOR_filltered2)
colnames(UMAP) = c("X","Y","chr","str","end","AB_signal_DAY0","AB_signal_DAY4","AB_signal_DIFF","DIFF_Q","Ets1","Ets1_Q","p300","p300_Q")

if(min_dist_para==0.0001){min_dist_para="00001"}
if(min_dist_para==0.001){min_dist_para="0001"}
if(min_dist_para==0.01){min_dist_para="001"}
if(min_dist_para==0.1){min_dist_para="01"}
if(min_dist_para==1){min_dist_para="1"}
if(min_dist_para==10){min_dist_para="10"}

if(spread_para==0.01){spread_para="001"}
if(spread_para==0.1){spread_para="01"}
if(spread_para==1){spread_para="1"}
if(spread_para==5){spread_para="5"}
if(spread_para==10){spread_para="10"}

path = paste("/Ets1/ProcessedData/UMAP",
             "_","DIFF",
             "_",metric_para,
             "_",n_neighbors_para,
             "_",n_components_para,
             "_",min_dist_para,
             "_",spread_para,
             ".txt",
             sep="")

write.table(UMAP, file = path, col.names=T, row.names=F, quote=F, sep="\t")
    } 

UMAP_viewer = function(inputPATH){

set.seed(7777)

UMAP_DIFF_correlation_50_2_1_1 = as.data.frame(vroom(inputPATH,delim = "\t",show_col_types = FALSE))
UMAP_DIFF_correlation_50_2_1_1_fit = UMAP_DIFF_correlation_50_2_1_1[,1:2]
UMAP_DIFF_correlation_50_2_1_1_label = paste(UMAP_DIFF_correlation_50_2_1_1[,3],
                                             UMAP_DIFF_correlation_50_2_1_1[,4],
                                             UMAP_DIFF_correlation_50_2_1_1[,5])
UMAP_DIFF_correlation_50_2_1_1_label = as.data.frame(UMAP_DIFF_correlation_50_2_1_1_label)
colnames(UMAP_DIFF_correlation_50_2_1_1_label) = "label"

result = kmeans(UMAP_DIFF_correlation_50_2_1_1_fit, centers = 22, iter.max = 10000, algorithm = "MacQueen")

UMAP_DIFF_correlation_50_2_1_1 = cbind(UMAP_DIFF_correlation_50_2_1_1,as.data.frame(result$cluster))
colnames(UMAP_DIFF_correlation_50_2_1_1)[ncol(UMAP_DIFF_correlation_50_2_1_1)] = "cluster"

options(repr.plot.width = 8, repr.plot.height = 8, repr.plot.res = 1000, repr.plot.pointsize = 30)

color_palette = c("darkblue", "#B71C1C","#283593", "#3949AB", "#5C6BC0", 
                  "#9FA8DA", "#757575", "#FFCDD2", "#E57373", "#E53935")

ggplot(UMAP_DIFF_correlation_50_2_1_1, aes(X, Y)) +
geom_point(aes(colour = factor(cluster)), size = 1) +
theme_classic(base_size = 12) +
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
    }

UMAP_DIFF_viewer = function(inputPATH){

UMAP = read.table(inputPATH,sep = "\t", header=T)

options(repr.plot.width = 8, repr.plot.height = 8, repr.plot.res = 1000, repr.plot.pointsize = 30)

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

UMAP_function(metric_para = "correlation",
              n_neighbors_para = 50,
              n_components_para = 2,
              min_dist_para = 1,
              spread_para = 1,
              n_threads_para = 30)

UMAP_viewer("/Ets1/ProcessedData/UMAP_DIFF_correlation_50_2_1_1.txt")

UMAP_DIFF_viewer("/Ets1/ProcessedData/UMAP_DIFF_correlation_50_2_1_1.txt")

