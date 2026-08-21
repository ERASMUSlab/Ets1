MATRIX_GM_UP = as.data.frame(vroom("/Ets1/ProcessedData/GM_MATRIX_UPPER.txt", delim = "\t",show_col_types = FALSE))
MATRIX_OM_UP = as.data.frame(vroom("/Ets1/ProcessedData/OM_MATRIX_UPPER.txt", delim = "\t",show_col_types = FALSE))

MATRIX_GM_UP_MAT = as.matrix(MATRIX_GM_UP)  
MATRIX_GM = MATRIX_GM_UP_MAT
MATRIX_GM[lower.tri(MATRIX_GM)] = t(MATRIX_GM)[lower.tri(MATRIX_GM)]
MATRIX_GM = as.data.frame(MATRIX_GM)

MATRIX_OM_UP_MAT = as.matrix(MATRIX_OM_UP)  
MATRIX_OM = MATRIX_OM_UP_MAT
MATRIX_OM[lower.tri(MATRIX_OM)] = t(MATRIX_OM)[lower.tri(MATRIX_OM)]
MATRIX_OM = as.data.frame(MATRIX_OM)

MATRIX_DIFF = MATRIX_OM - MATRIX_GM

options(repr.plot.width = 5, repr.plot.height = 5, repr.plot.res = 1000, repr.plot.pointsize = 40)

MAP = MATRIX_DIFF

conte=MAP
rownames(conte)=c(1:nrow(MAP))

PCA=conte
fig = pheatmap(PCA,
               cluster_rows = F,
               cluster_cols = F,
               show_rownames = F,
               show_colnames = F,
               legend = FALSE,
               border_color = FALSE,
               color = colorRampPalette(c("darkblue","white","darkred"))(100),
               breaks=seq(from=-3,to=3,length.out = 100))

