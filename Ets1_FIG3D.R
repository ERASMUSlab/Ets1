TOP_BOTTOM_MAT = function(chr_para,ratio,outpath,histone_input){

chr_DF = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
           "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19")

len_DF = c(195471971,182113224,160039680,156508116,151834684,149736546,145441459,129401213,124595110,130694993,
           122082543,120129022,120421639,124902244,104043685,98207768,94987271,90702639,61431566)

para_DF = cbind(chr_DF,len_DF)
para_DF = as.data.frame(para_DF)

len = as.numeric(para_DF[which(para_DF[,1] == chr_para),2])

DAY0_WT_H3K27ac_region = read.table(histone_input,sep = "\t", header=F,fill = TRUE)
DAY0_WT_H3K27ac_region = DAY0_WT_H3K27ac_region[,1:3]
colnames(DAY0_WT_H3K27ac_region) = c("chr","starts","ends")

DAY0_WT_H3K27ac_region = DAY0_WT_H3K27ac_region[which(DAY0_WT_H3K27ac_region[,1] == chr_para),]

ABcompartmentDATA = read.table("/data3/psg/Ets1/ProcessedData/ABcompartmentDATA.txt",sep = "\t", header=F,fill = TRUE)
ABcompartmentDATA$label = paste(ABcompartmentDATA[,1],ABcompartmentDATA[,2],ABcompartmentDATA[,3])
ABcompartmentDATA = ABcompartmentDATA[,4:5]

ABcompartment_100kbDATA = read.table("/data3/psg/Ets1/ProcessedData/ABcompartment_100kbDATA.txt",sep = "\t", header=T,fill = TRUE)
ABcompartment_100kbDATA = inner_join(ABcompartment_100kbDATA,ABcompartmentDATA,by="label")

Acompartment = ABcompartment_100kbDATA[which(ABcompartment_100kbDATA[,8] == "GM_A"),1:3]
Bcompartment = ABcompartment_100kbDATA[which(ABcompartment_100kbDATA[,8] == "GM_B"),1:3]

Acompartment = Acompartment[which(Acompartment[,1] == chr_para),]
Bcompartment = Bcompartment[which(Bcompartment[,1] == chr_para),]

Acompartment[,2] = Acompartment[,2] -1
Bcompartment[,2] = Bcompartment[,2] -1

Acompartment = Acompartment[order(Acompartment[,2]),]
Bcompartment = Bcompartment[order(Bcompartment[,2]),]

Acompartment[,2] = as.integer(Acompartment[,2])
Acompartment[,3] = as.integer(Acompartment[,3])

Bcompartment[,2] = as.integer(Bcompartment[,2])
Bcompartment[,3] = as.integer(Bcompartment[,3])

Acompartment$Asite_label = paste(Acompartment[,1],
                                      Acompartment[,2],
                                      Acompartment[,3])

Acompartment$Bsite_label = paste(Acompartment[,1],
                                      Acompartment[,2],
                                      Acompartment[,3])

Bcompartment$Asite_label = paste(Bcompartment[,1],
                                      Bcompartment[,2],
                                      Bcompartment[,3])

Bcompartment$Bsite_label = paste(Bcompartment[,1],
                                      Bcompartment[,2],
                                      Bcompartment[,3])

start = seq(0,len,by=100000)
start = as.integer(start)
start = as.data.frame(start)
chrname = c(1:nrow(start))

DAY0_WT_corrected_binlist_FULL = cbind(chrname,start)
DAY0_WT_corrected_binlist_FULL[,1] = chr_para
colnames(DAY0_WT_corrected_binlist_FULL) = c("chrA","startA")
DAY0_WT_corrected_binlist_FULL$endA = DAY0_WT_corrected_binlist_FULL[,2]+100000
DAY0_WT_corrected_binlist_FULL[,3] = as.integer(DAY0_WT_corrected_binlist_FULL[,3])

DAY0_WT_corrected_binlist_DB = DAY0_WT_corrected_binlist_FULL
DAY0_WT_corrected_binlist_DB$chrB = 0
DAY0_WT_corrected_binlist_DB$startB = 0
DAY0_WT_corrected_binlist_DB$endB = 0

DAY0_WT_corrected_binlist_DB$chrB = DAY0_WT_corrected_binlist_FULL[1,1]
DAY0_WT_corrected_binlist_DB$startB = DAY0_WT_corrected_binlist_FULL[1,2]
DAY0_WT_corrected_binlist_DB$endB = DAY0_WT_corrected_binlist_FULL[1,3]

for(i in 2:nrow(DAY0_WT_corrected_binlist_FULL)){
    DAY0_WT_corrected_binlist_DB_frag = DAY0_WT_corrected_binlist_FULL
    DAY0_WT_corrected_binlist_DB_frag$chrB = 0
    DAY0_WT_corrected_binlist_DB_frag$startB = 0
    DAY0_WT_corrected_binlist_DB_frag$endB = 0
    
    DAY0_WT_corrected_binlist_DB_frag$chrB = DAY0_WT_corrected_binlist_FULL[i,1]
    DAY0_WT_corrected_binlist_DB_frag$startB = DAY0_WT_corrected_binlist_FULL[i,2]
    DAY0_WT_corrected_binlist_DB_frag$endB = DAY0_WT_corrected_binlist_FULL[i,3]

    DAY0_WT_corrected_binlist_DB = rbind(DAY0_WT_corrected_binlist_DB,DAY0_WT_corrected_binlist_DB_frag)
    }

DAY0_WT_corrected_binAC = DAY0_WT_corrected_binlist_FULL
DAY0_WT_corrected_binAC[,4] = 0

for(i in 1:nrow(DAY0_WT_corrected_binAC)){
    DAY0_WT_H3K27ac_region_tmp1 = DAY0_WT_H3K27ac_region[which(DAY0_WT_H3K27ac_region[,2] < DAY0_WT_corrected_binAC[i,3]),]
    DAY0_WT_H3K27ac_region_type1 = DAY0_WT_H3K27ac_region_tmp1[which(DAY0_WT_H3K27ac_region_tmp1[,3] > DAY0_WT_corrected_binAC[i,2]),]
    
    overlaped = DAY0_WT_H3K27ac_region_type1
    overlaped = unique(overlaped)

    if(nrow(overlaped)>0){
        DAY0_WT_corrected_binAC[i,4] = nrow(overlaped)
        }
}

colnames(DAY0_WT_corrected_binAC)[4] = "acA"

DAY0_WT_corrected_binAC_AB = DAY0_WT_corrected_binAC
DAY0_WT_corrected_binAC_AB$chrB = 0
DAY0_WT_corrected_binAC_AB$startB = 0
DAY0_WT_corrected_binAC_AB$endB = 0
DAY0_WT_corrected_binAC_AB$acB = 0 

DAY0_WT_corrected_binAC_AB$chrB = DAY0_WT_corrected_binAC[1,1]
DAY0_WT_corrected_binAC_AB$startB = DAY0_WT_corrected_binAC[1,2]
DAY0_WT_corrected_binAC_AB$endB = DAY0_WT_corrected_binAC[1,3]
DAY0_WT_corrected_binAC_AB$acB = DAY0_WT_corrected_binAC[1,4]

for(i in 2:nrow(DAY0_WT_corrected_binAC)){
    DAY0_WT_corrected_binAC_AB_frag = DAY0_WT_corrected_binAC
    DAY0_WT_corrected_binAC_AB_frag$chrB = 0
    DAY0_WT_corrected_binAC_AB_frag$startB = 0
    DAY0_WT_corrected_binAC_AB_frag$endB = 0
    DAY0_WT_corrected_binAC_AB_frag$acB = 0 
    
    DAY0_WT_corrected_binAC_AB_frag$chrB = DAY0_WT_corrected_binAC[i,1]
    DAY0_WT_corrected_binAC_AB_frag$startB = DAY0_WT_corrected_binAC[i,2]
    DAY0_WT_corrected_binAC_AB_frag$endB = DAY0_WT_corrected_binAC[i,3]
    DAY0_WT_corrected_binAC_AB_frag$acB = DAY0_WT_corrected_binAC[i,4]

    DAY0_WT_corrected_binAC_AB = rbind(DAY0_WT_corrected_binAC_AB,DAY0_WT_corrected_binAC_AB_frag)
    }

DAY0_bin_pAC = DAY0_WT_corrected_binAC_AB
DAY0_bin_pAC$pAC = DAY0_bin_pAC[,4]+DAY0_bin_pAC[,8]

DAY0_bin_pAC_topONE_tmp = DAY0_bin_pAC[order(-DAY0_bin_pAC[,9]),]
DAY0_bin_pAC_topONE_tmp[round(nrow(DAY0_bin_pAC)*ratio):nrow(DAY0_bin_pAC_topONE_tmp),9] = 0

DAY0_bin_pAC_topONE_tmp$label = paste(DAY0_bin_pAC_topONE_tmp[,1],
                                            DAY0_bin_pAC_topONE_tmp[,2],
                                            DAY0_bin_pAC_topONE_tmp[,3],
                                            DAY0_bin_pAC_topONE_tmp[,5],
                                            DAY0_bin_pAC_topONE_tmp[,6],
                                            DAY0_bin_pAC_topONE_tmp[,7])

DAY0_bin_pAC_topONE = DAY0_bin_pAC
DAY0_bin_pAC_topONE$label = paste(DAY0_bin_pAC_topONE[,1],
                                        DAY0_bin_pAC_topONE[,2],
                                        DAY0_bin_pAC_topONE[,3],
                                        DAY0_bin_pAC_topONE[,5],
                                        DAY0_bin_pAC_topONE[,6],
                                        DAY0_bin_pAC_topONE[,7])

DAY0_bin_pAC_topONE_tmp = DAY0_bin_pAC_topONE_tmp[,c(10,9)]
DAY0_bin_pAC_topONE = DAY0_bin_pAC_topONE[,-9]
DAY0_bin_pAC_topONE = inner_join(DAY0_bin_pAC_topONE,DAY0_bin_pAC_topONE_tmp,by="label")

DAY0_bin_pAC_bottomONE_tmp = DAY0_bin_pAC[order(-DAY0_bin_pAC[,9]),]
DAY0_bin_pAC_bottomONE_tmp[1:(nrow(DAY0_bin_pAC_bottomONE_tmp)-round(nrow(DAY0_bin_pAC)*ratio)),9] = 0

DAY0_bin_pAC_bottomONE_tmp$label = paste(DAY0_bin_pAC_bottomONE_tmp[,1],
                                            DAY0_bin_pAC_bottomONE_tmp[,2],
                                            DAY0_bin_pAC_bottomONE_tmp[,3],
                                            DAY0_bin_pAC_bottomONE_tmp[,5],
                                            DAY0_bin_pAC_bottomONE_tmp[,6],
                                            DAY0_bin_pAC_bottomONE_tmp[,7])

DAY0_bin_pAC_bottomONE = DAY0_bin_pAC
DAY0_bin_pAC_bottomONE$label = paste(DAY0_bin_pAC_bottomONE[,1],
                                        DAY0_bin_pAC_bottomONE[,2],
                                        DAY0_bin_pAC_bottomONE[,3],
                                        DAY0_bin_pAC_bottomONE[,5],
                                        DAY0_bin_pAC_bottomONE[,6],
                                        DAY0_bin_pAC_bottomONE[,7])

DAY0_bin_pAC_bottomONE_tmp = DAY0_bin_pAC_bottomONE_tmp[,c(10,9)]
DAY0_bin_pAC_bottomONE = DAY0_bin_pAC_bottomONE[,-9]
DAY0_bin_pAC_bottomONE = inner_join(DAY0_bin_pAC_bottomONE,DAY0_bin_pAC_bottomONE_tmp,by="label")
DAY0_bin_pAC_bottomONE[,10] = DAY0_bin_pAC_bottomONE[,10]*(-1)

rowBIN = seq(1:nrow(start))
BINnum = rep(rowBIN, times = nrow(start))
BINnum = as.data.frame(BINnum)

DAY0_bin_pAC_BINnum = cbind(DAY0_bin_pAC,BINnum)

DAY0_bin_pAC_MAT = DAY0_bin_pAC_BINnum[which(DAY0_bin_pAC_BINnum$BINnum == 1),8]
DAY0_bin_pAC_MAT = as.data.frame(DAY0_bin_pAC_MAT)
colnames(DAY0_bin_pAC_MAT) = NULL

for(i in 2:nrow(start)){
    DAY0_bin_pAC_MATfrag = DAY0_bin_pAC_BINnum[which(DAY0_bin_pAC_BINnum$BINnum == i),8]
    DAY0_bin_pAC_MATfrag = as.data.frame(DAY0_bin_pAC_MATfrag)
    colnames(DAY0_bin_pAC_MATfrag) = NULL
    DAY0_bin_pAC_MAT = cbind(DAY0_bin_pAC_MAT,DAY0_bin_pAC_MATfrag)
    }

colnames(DAY0_bin_pAC_MAT) = NULL

DAY0_bin_pAC_MAT = as.matrix(DAY0_bin_pAC_MAT)
pmean <- function(x,y) (x+y)/2
DAY0_bin_pAC_MAT[] <- pmean(DAY0_bin_pAC_MAT, matrix(DAY0_bin_pAC_MAT, nrow(DAY0_bin_pAC_MAT), byrow=TRUE))     

rowBIN = seq(1:nrow(start))
BINnum = rep(rowBIN, times = nrow(start))
BINnum = as.data.frame(BINnum)

DAY0_bin_pAC_topONE_BINnum = cbind(DAY0_bin_pAC_topONE,BINnum)

DAY0_bin_pAC_topONE_MAT = DAY0_bin_pAC_topONE_BINnum[which(DAY0_bin_pAC_topONE_BINnum$BINnum == 1),10]
DAY0_bin_pAC_topONE_MAT = as.data.frame(DAY0_bin_pAC_topONE_MAT)
colnames(DAY0_bin_pAC_topONE_MAT) = NULL

for(i in 2:nrow(start)){
    DAY0_bin_pAC_topONE_MATfrag = DAY0_bin_pAC_topONE_BINnum[which(DAY0_bin_pAC_topONE_BINnum$BINnum == i),10]
    DAY0_bin_pAC_topONE_MATfrag = as.data.frame(DAY0_bin_pAC_topONE_MATfrag)
    colnames(DAY0_bin_pAC_topONE_MATfrag) = NULL
    DAY0_bin_pAC_topONE_MAT = cbind(DAY0_bin_pAC_topONE_MAT,DAY0_bin_pAC_topONE_MATfrag)
    }

colnames(DAY0_bin_pAC_topONE_MAT) = NULL

DAY0_bin_pAC_topONE_MAT = as.matrix(DAY0_bin_pAC_topONE_MAT)
pmean <- function(x,y) (x+y)/2
DAY0_bin_pAC_topONE_MAT[] <- pmean(DAY0_bin_pAC_topONE_MAT, matrix(DAY0_bin_pAC_topONE_MAT, nrow(DAY0_bin_pAC_topONE_MAT), byrow=TRUE))

rowBIN = seq(1:nrow(start))
BINnum = rep(rowBIN, times = nrow(start))
BINnum = as.data.frame(BINnum)

DAY0_bin_pAC_bottomONE_BINnum = cbind(DAY0_bin_pAC_bottomONE,BINnum)

DAY0_bin_pAC_bottomONE_MAT = DAY0_bin_pAC_bottomONE_BINnum[which(DAY0_bin_pAC_bottomONE_BINnum$BINnum == 1),10]
DAY0_bin_pAC_bottomONE_MAT = as.data.frame(DAY0_bin_pAC_bottomONE_MAT)
colnames(DAY0_bin_pAC_bottomONE_MAT) = NULL

for(i in 2:nrow(start)){
    DAY0_bin_pAC_bottomONE_MATfrag = DAY0_bin_pAC_bottomONE_BINnum[which(DAY0_bin_pAC_bottomONE_BINnum$BINnum == i),10]
    DAY0_bin_pAC_bottomONE_MATfrag = as.data.frame(DAY0_bin_pAC_bottomONE_MATfrag)
    colnames(DAY0_bin_pAC_bottomONE_MATfrag) = NULL
    DAY0_bin_pAC_bottomONE_MAT = cbind(DAY0_bin_pAC_bottomONE_MAT,DAY0_bin_pAC_bottomONE_MATfrag)
    }

colnames(DAY0_bin_pAC_bottomONE_MAT) = NULL

DAY0_bin_pAC_bottomONE_MAT = as.matrix(DAY0_bin_pAC_bottomONE_MAT)
pmean <- function(x,y) (x+y)/2
DAY0_bin_pAC_bottomONE_MAT[] <- pmean(DAY0_bin_pAC_bottomONE_MAT, matrix(DAY0_bin_pAC_bottomONE_MAT, nrow(DAY0_bin_pAC_bottomONE_MAT), byrow=TRUE)) 


DAY0_bin_pAC_TEN_MAT = DAY0_bin_pAC_topONE_MAT+DAY0_bin_pAC_bottomONE_MAT

DAY0_bin_pAC_bottomTEN = DAY0_bin_pAC_bottomONE_BINnum[which(DAY0_bin_pAC_bottomONE_BINnum$pAC < 0),c(1,2,3,5,6,7)]

DAY0_bin_pAC_bottomTEN$Asite_label = paste(DAY0_bin_pAC_bottomTEN[,1],
                                                DAY0_bin_pAC_bottomTEN[,2],
                                                DAY0_bin_pAC_bottomTEN[,3])

DAY0_bin_pAC_bottomTEN$Bsite_label = paste(DAY0_bin_pAC_bottomTEN[,4],
                                                DAY0_bin_pAC_bottomTEN[,5],
                                                DAY0_bin_pAC_bottomTEN[,6])

DAY0_bin_pAC_bottomTEN_Asite = DAY0_bin_pAC_bottomTEN

DAY0_bin_pAC_bottomTEN_Asite$Asite = "NONE"

for(i in 1:nrow(DAY0_bin_pAC_bottomTEN_Asite)){
    Acompartment_corr = Acompartment[which(Acompartment[,4] == DAY0_bin_pAC_bottomTEN_Asite[i,7]),]
    Bcompartment_corr = Bcompartment[which(Bcompartment[,4] == DAY0_bin_pAC_bottomTEN_Asite[i,7]),]

    if(nrow(Acompartment_corr)>0){
        DAY0_bin_pAC_bottomTEN_Asite[i,9] = "A"
        }

    if(nrow(Bcompartment_corr)>0){
        DAY0_bin_pAC_bottomTEN_Asite[i,9] = "B"
        }
}

DAY0_bin_pAC_bottomTEN_FULL = DAY0_bin_pAC_bottomTEN_Asite

DAY0_bin_pAC_bottomTEN_FULL$Bsite = "NONE"

for(i in 1:nrow(DAY0_bin_pAC_bottomTEN_FULL)){
    Acompartment_corr = Acompartment[which(Acompartment[,4] == DAY0_bin_pAC_bottomTEN_FULL[i,8]),]
    Bcompartment_corr = Bcompartment[which(Bcompartment[,4] == DAY0_bin_pAC_bottomTEN_FULL[i,8]),]

    if(nrow(Acompartment_corr)>0){
        DAY0_bin_pAC_bottomTEN_FULL[i,10] = "A"
        }

    if(nrow(Bcompartment_corr)>0){
        DAY0_bin_pAC_bottomTEN_FULL[i,10] = "B"
        }
}

DAY0_bin_pAC_topTEN = DAY0_bin_pAC_topONE_BINnum[which(DAY0_bin_pAC_topONE_BINnum$pAC > 0),c(1,2,3,5,6,7)]

DAY0_bin_pAC_topTEN$Asite_label = paste(DAY0_bin_pAC_topTEN[,1],
                                                DAY0_bin_pAC_topTEN[,2],
                                                DAY0_bin_pAC_topTEN[,3])

DAY0_bin_pAC_topTEN$Bsite_label = paste(DAY0_bin_pAC_topTEN[,4],
                                                DAY0_bin_pAC_topTEN[,5],
                                                DAY0_bin_pAC_topTEN[,6])

DAY0_bin_pAC_topTEN_Asite = DAY0_bin_pAC_topTEN

DAY0_bin_pAC_topTEN_Asite$Asite = "NONE"

for(i in 1:nrow(DAY0_bin_pAC_topTEN_Asite)){
    Acompartment_corr = Acompartment[which(Acompartment[,4] == DAY0_bin_pAC_topTEN_Asite[i,7]),]
    Bcompartment_corr = Bcompartment[which(Bcompartment[,4] == DAY0_bin_pAC_topTEN_Asite[i,7]),]

    if(nrow(Acompartment_corr)>0){
        DAY0_bin_pAC_topTEN_Asite[i,9] = "A"
        }

    if(nrow(Bcompartment_corr)>0){
        DAY0_bin_pAC_topTEN_Asite[i,9] = "B"
        }
}

DAY0_bin_pAC_topTEN_FULL = DAY0_bin_pAC_topTEN_Asite

DAY0_bin_pAC_topTEN_FULL$Bsite = "NONE"

for(i in 1:nrow(DAY0_bin_pAC_topTEN_FULL)){
    Acompartment_corr = Acompartment[which(Acompartment[,4] == DAY0_bin_pAC_topTEN_FULL[i,8]),]
    Bcompartment_corr = Bcompartment[which(Bcompartment[,4] == DAY0_bin_pAC_topTEN_FULL[i,8]),]

    if(nrow(Acompartment_corr)>0){
        DAY0_bin_pAC_topTEN_FULL[i,10] = "A"
        }

    if(nrow(Bcompartment_corr)>0){
        DAY0_bin_pAC_topTEN_FULL[i,10] = "B"
        }
}
    
write.table(DAY0_bin_pAC_TEN_MAT, file = outpath ,col.names=T, row.names=F, quote=F, sep="\t")
    }

TOP_BOTTOM_MAT(chr_para = "chr5", ratio = 0.1, 
               histone_input="/data3/psg/Ets1/ProcessedData/GM_Histone_fit.bedpe",
               outpath = "/data3/psg/Ets1/ProcessedData/Histone_chr5_TOPBOTTOM_10_MAP.txt")

TOP_BOTTOM_MAT(chr_para = "chr5", ratio = 0.1, 
               histone_input="/data3/psg/Ets1/ProcessedData/GM_H3K27ac_fit.bedpe",
               outpath = "/data3/psg/Ets1/ProcessedData/H3K27ac_chr5_TOPBOTTOM_10_MAP.txt")

TOP_BOTTOM_MAT(chr_para = "chr5", ratio = 0.1, 
               histone_input="/data3/psg/Ets1/ProcessedData/GM_H3K4me1_fit.bedpe",
               outpath = "/data3/psg/Ets1/ProcessedData/H3K4me1_chr5_TOPBOTTOM_10_MAP.txt")

TOP_BOTTOM_MAT(chr_para = "chr5", ratio = 0.1, 
               histone_input="/data3/psg/Ets1/ProcessedData/GM_H3K4me3_fit.bedpe",
               outpath = "/data3/psg/Ets1/ProcessedData/H3K4me3_chr5_TOPBOTTOM_10_MAP.txt")

dataINPUT = read.table("/data3/psg/Ets1/ProcessedData/Histone_chr5_TOPBOTTOM_10_MAP.txt",sep = "\t", header=T,fill = TRUE)

options(repr.plot.width = 5, repr.plot.height = 5, repr.plot.res = 1000, repr.plot.pointsize = 10)

data = dataINPUT

conte=data
namen = c(1:nrow(data))
rownames(conte)=as.vector(namen)
PCA=conte
pheatmap(PCA,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = F,
         legend = FALSE,
         border_color = FALSE,
         color = colorRampPalette(c("darkblue","white","darkred"))(100),
         breaks=seq(from=-1,to=1,length.out = 100))

H3K27ac_chr5_MAP = read.table("/data3/psg/Ets1/ProcessedData/H3K27ac_chr5_MAP.txt",sep = "\t", header=T,fill = TRUE)
H3K4me1_chr5_MAP = read.table("/data3/psg/Ets1/ProcessedData/H3K4me1_chr5_MAP.txt",sep = "\t", header=T,fill = TRUE)
H3K4me3_chr5_MAP = read.table("/data3/psg/Ets1/ProcessedData/H3K4me3_chr5_MAP.txt",sep = "\t", header=T,fill = TRUE)

chr5_MAP = (H3K27ac_chr5_MAP+H3K4me1_chr5_MAP+H3K4me3_chr5_MAP)/3

data = chr5_MAP

options(repr.plot.width = 5, repr.plot.height = 5, repr.plot.res = 1000, repr.plot.pointsize = 10)
conte=data
namen = c(1:nrow(data))
rownames(conte)=as.vector(namen)
PCA=conte
pheatmap(PCA,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = F,
         legend = FALSE,
         border_color = FALSE,
         color = colorRampPalette(c("white","white","tan1","tan1","tan2","tan3","tan4","tan4"))(100),
         breaks=seq(from=0,to=1600,length.out = 100))

