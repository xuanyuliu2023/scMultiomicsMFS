library(Matrix)
library(Seurat)
library(ggplot2)
STdata <- get(load('STobj.GEXintegrated.Rdata'))
DefaultAssay(STdata) <- 'SCT'
feature='FOXN3'
SpatialFeaturePlot(STdata, slot = 'data',features = feature, pt.size.factor = 1.6,alpha = c(0.1, 1) ) + theme(legend.position = "right")
ggsave(file=paste0(feature,'_spatialFearturePlot.alpha0.1.jpg'),dpi = 600,width=12.5,height=9.06,units="in")
ggsave(file=paste0(feature,'_spatialFearturePlot.alpha0.1.pdf'))


#SCSE prep
normData <- t(as.matrix(STdata@assays$SCT@data))
write.table( normData, "Normalized.expr.tsv", sep="\t", row.names = TRUE, col.names=NA)


df <- read.table(file='c2_CP_HUMAN_v7.2_Normalized.expr.tsv',header=T,check.names =F,sep='\t')
df2 <- df[,2:ncol(df)]
row.names(df2) <- df$id
pathwayActivityMatrix <- as(t(df2), "dgCMatrix")
STdata@assays$SCT@data <- pathwayActivityMatrix

feature='WP_PI3KAKT_SIGNALING_PATHWAY'
feature='KEGG_TGF_BETA_SIGNALING_PATHWAY'
feature='ST_ERK1_ERK2_MAPK_PATHWAY'
feature='REACTOME_ERK_MAPK_TARGETS'
#pt.size.factor- This will scale the size of the spots. Default is 1.6
#alpha - minimum and maximum transparency. Default is c(1, 1).
SpatialFeaturePlot(STdata, features = feature, slot = 'data', pt.size.factor = 1.6,alpha = c(1, 1),stroke=0.1 ) + theme(legend.position = "top")
ggsave(file=paste0(feature,'_spatialFearturePlot.alpha1.pdf'))

ggsave(file=paste0(feature,'_spatialFearturePlot.alpha1.jpg'),dpi = 600,width=17.4,height=6.71,units="in")


