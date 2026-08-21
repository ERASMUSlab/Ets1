compartment_heatmap = function(input){
    
options(repr.plot.width = 1.5, repr.plot.height = 2.5, repr.plot.res = 300, repr.plot.pointsize = 3)

conte=input[,4:5]
rownames(conte)=input[,7]
PCA=conte
pheatmap(PCA,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = F,
         colorRampPalette(c("yellow3", "black", "dodgerblue3"))(100),
         legend = FALSE,
         border_color =NA,
         breaks=seq(from=-0.05,to=0.05,length.out = 100))
    } 

total_100kb_bin = read.table("/Ets1/ProcessedData/ABcompartment_100kbDATA.txt",sep = "\t", header=T,fill = TRUE)

A_100kb_bin = total_100kb_bin[which(total_100kb_bin$DAY0 > 0),]

strA_100kb_bin = A_100kb_bin[which(A_100kb_bin$DAY0 > 0.00654203),]
strA_100kb_bin = strA_100kb_bin[which(strA_100kb_bin$DIFF > 0.002984229),]

stableA_100kb_bin = A_100kb_bin[which(A_100kb_bin$DAY0 > 0.00654203),]
stableA_100kb_bin = stableA_100kb_bin[which(stableA_100kb_bin$DIFF <= 0.002984229),]
stableA_100kb_bin = stableA_100kb_bin[which(stableA_100kb_bin$DAY4 >= 0.00654203),]

B_100kb_bin = total_100kb_bin[which(total_100kb_bin$DAY0 < 0),]

strB_100kb_bin = B_100kb_bin[which(B_100kb_bin$DAY0 < -0.006892938),]
strB_100kb_bin = strB_100kb_bin[which(strB_100kb_bin$DIFF < -0.003471521),]

stableB_100kb_bin = B_100kb_bin[which(B_100kb_bin$DAY0 < -0.006892938),]
stableB_100kb_bin = stableB_100kb_bin[which(stableB_100kb_bin$DIFF >= -0.003471521),]
stableB_100kb_bin = stableB_100kb_bin[which(stableB_100kb_bin$DAY4 <= -0.006892938),]

changeA_100kb_bin = total_100kb_bin[which(total_100kb_bin$DAY0 > 0.00654203),]
changeA_100kb_bin = changeA_100kb_bin[which(changeA_100kb_bin$DAY4 < 0),]

changeB_100kb_bin = total_100kb_bin[which(total_100kb_bin$DAY0 < -0.006892938),]
changeB_100kb_bin = changeB_100kb_bin[which(changeB_100kb_bin$DAY4 > 0),]


subtle_100kb_bin_tmp = rbind(strA_100kb_bin,
                             stableA_100kb_bin,
                             changeA_100kb_bin,
                             strB_100kb_bin,
                             stableB_100kb_bin,
                             changeB_100kb_bin)
subtle_100kb_bin_tmp = as.data.frame(subtle_100kb_bin_tmp[,7])
colnames(subtle_100kb_bin_tmp) = "label"
subtle_100kb_bin_tmp = unique(subtle_100kb_bin_tmp)

total_100kb_bin_label = as.data.frame(total_100kb_bin[,7])
colnames(total_100kb_bin_label) = "label"
total_100kb_bin_label = unique(total_100kb_bin_label)

subtle_100kb_bin_tmp = anti_join(total_100kb_bin_label,subtle_100kb_bin_tmp,by="label")
subtle_100kb_bin = inner_join(total_100kb_bin,subtle_100kb_bin_tmp,by="label")

round(nrow(strA_100kb_bin)/nrow(total_100kb_bin)*100,2)
compartment_heatmap(input = strA_100kb_bin)

round(nrow(strB_100kb_bin)/nrow(total_100kb_bin)*100,2)
compartment_heatmap(input = strB_100kb_bin)

round(nrow(stableA_100kb_bin)/nrow(total_100kb_bin)*100,2)
compartment_heatmap(input = stableA_100kb_bin)

round(nrow(stableB_100kb_bin)/nrow(total_100kb_bin)*100,2)
compartment_heatmap(input = stableB_100kb_bin)

round(nrow(changeA_100kb_bin)/nrow(total_100kb_bin)*100,2)
compartment_heatmap(input = changeA_100kb_bin)

round(nrow(changeB_100kb_bin)/nrow(total_100kb_bin)*100,2)
compartment_heatmap(input = changeB_100kb_bin)

round(nrow(subtle_100kb_bin)/nrow(total_100kb_bin)*100,2)
compartment_heatmap(input = subtle_100kb_bin)

