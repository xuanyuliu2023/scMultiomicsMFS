load('/work/project/xuanyu/project/MFSmultiome/Signac/Integration/iter5/step2/seuObj.step2.Rdata')
library(Signac)
library(Seurat)
library(future)
plan("multicore", workers = 5)
options(future.globals.maxSize = 200 * 1024 ^ 3)
plan()
library(plyr)
MotifNames=names(seuObj@assays$peaks@motifs@motif.names)
TFs=as.character(seuObj@assays$peaks@motifs@motif.names)

DefaultAssay(seuObj) <- 'chromvar'
seuObj@meta.data$cellTypeGroup <- paste0(seuObj@meta.data$cellType,"_",seuObj@meta.data$group)
Idents(seuObj) <- 'cellTypeGroup'

DiffMotif_test <- function(targetCluster){
diffMotifs_diseased_vs_healthy <- FindMarkers(
  object = seuObj,
  ident.1 = paste0(targetCluster,'_diseased'),
  ident.2 = paste0(targetCluster,'_healthy'),
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)
diffMotifs_diseased_vs_healthy$peak <-rownames(diffMotifs_diseased_vs_healthy)
diffMotifs_diseased_vs_healthy$TF <- mapvalues(diffMotifs_diseased_vs_healthy$peak,from=MotifNames, to=TFs)
diffMotifs_diseased_vs_healthy <- subset(diffMotifs_diseased_vs_healthy,p_val_adj <0.05)

write.table(diffMotifs_diseased_vs_healthy,file=paste0('diffMotifs_',targetCluster,'_diseased_vs_',targetCluster,'_healthy.tsv'),sep='\t',col.names=T,row.names=F,quote=F)
}

cellTypeGroups <- unique(seuObj@meta.data$cellType)
lapply(cellTypeGroups,DiffMotif_test)

FeaturePlot(
  object = seuObj,
  features = "MA0113.3",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)

MotifPlot(
  object = seuObj,
  motifs = "MA0113.3",
  assay = 'peaks'
)
