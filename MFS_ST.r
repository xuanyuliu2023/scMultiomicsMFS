library(Signac)
library(Seurat)
library(future)
library(plyr)
plan("multicore", workers = 10)
options(future.globals.maxSize = 200 * 1024 ^ 3)
plan()

data_dir <- '/work/project/xuanyu/project/MFSmultiome/ST/spaceranger/AD2242/outs/'
STdata.AD2242 <- Load10X_Spatial(data.dir = data_dir)
STdata.AD2242@meta.data$sampleID <-rep('AD2242',ncol(STdata.AD2242))

data_dir <- '/work/project/xuanyu/project/MFSmultiome/ST/spaceranger/AD2287/outs/'
STdata.AD2287 <- Load10X_Spatial(data.dir = data_dir)
STdata.AD2287@meta.data$sampleID <-rep('AD2287',ncol(STdata.AD2287))

data_dir <- '/work/project/xuanyu/project/MFSmultiome/ST/spaceranger/YZ75/outs/'
STdata.YZ75 <- Load10X_Spatial(data.dir = data_dir)
STdata.YZ75@meta.data$sampleID <-rep('YZ75',ncol(STdata.YZ75))

data_dir <- '/work/project/xuanyu/project/MFSmultiome/ST/spaceranger/YZ74/outs/'
STdata.YZ74 <- Load10X_Spatial(data.dir = data_dir)
STdata.YZ74@meta.data$sampleID <-rep('YZ74',ncol(STdata.YZ74))

seuObj.list <- c(AD2242=STdata.AD2242,AD2287=STdata.AD2287,YZ75=STdata.YZ75,YZ74=STdata.YZ74)
# ---------- Gene expression data processing ------------------
# Normalize datasets individually by SCTransform(), instead of NormalizeData() prior to integration
# SCTransform (Hafemeister and Satija, 2019) procedure omits the need for heuristic steps including pseudocount addition or log-transformation
# The latest version of sctransform also supports using glmGamPoi package which substantially improves the speed of the learning procedure
#ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
# SCTransform is designed to eliminate the influence of sequencing depth. The nCount_RNA is used to construct the sct mode, so you don't need to put it into vars.to.regress.
seuObj.list <- lapply(X = seuObj.list, function(i){SCTransform(i, method = "glmGamPoi",assay='Spatial')})
# Typically use 3,000 or more features for analysis downstream of sctransform
features <- SelectIntegrationFeatures(object.list = seuObj.list, nfeatures = 3000, assay =c('SCT','SCT','SCT','SCT'))
# Run the PrepSCTIntegration() function prior to identifying anchors
seuObj.list <- PrepSCTIntegration(object.list = seuObj.list, anchor.features = features, assay =c('SCT','SCT','SCT','SCT'))

# integration
GEX.anchors <- FindIntegrationAnchors(object.list = seuObj.list, normalization.method = "SCT", anchor.features = features)
STobj.GEXintegrated <- IntegrateData(anchorset = GEX.anchors, normalization.method = "SCT",new.assay.name = "GEXintegrated",dims = 1:20)

#Perform dimensionality reduction by PCA and UMAP embedding
STobj.GEXintegrated <- RunPCA(STobj.GEXintegrated, verbose = FALSE, assay = "GEXintegrated",npcs = 50)
STobj.GEXintegrated <- RunUMAP(STobj.GEXintegrated, reduction = "pca", assay = "GEXintegrated", dims = 1:20)
STobj.GEXintegrated <- FindNeighbors(STobj.GEXintegrated, reduction = "pca", dims = 1:20)
STobj.GEXintegrated <- FindClusters(STobj.GEXintegrated, resolution = 0.6)
STobj.GEXintegrated <- FindClusters(STobj.GEXintegrated, resolution = 0.4)
STobj.GEXintegrated <- FindClusters(STobj.GEXintegrated, resolution = 0.3)
STobj.GEXintegrated <- FindClusters(STobj.GEXintegrated, resolution = 0.2)

STobj.GEXintegrated$group <- mapvalues(STobj.GEXintegrated$sampleID,from=c('AD2242','AD2287','YZ75','YZ74'),
to=c('diseased','diseased','healthy','healthy'))
# visualization
# Dimplot
library(ggplot2)

# UMAP categary data plot

categary.color.pallet  <- as.character(c('#1f77b4','#ff7f0e','#2ca02c','#d62728','#7733B7',
  '#8c564b','#E236AF','#8C8C00','#700B35','#06BACE',
"#193C3E", "#992756", "#DADD00", "#F76AA7", "#C4796A",
  "#246A73", "#FF6360", "#4DDD30", "#FFA449", "#74A6E8",
  "#955BC1","#FF0044","#F6757A","#265C42","#00D8FF",
  "#63C74D","#800026","#bd0026"))
categary='GEXintegrated_snn_res.0.4'
categary='sampleID'
categary='group'
DimPlot(STobj.GEXintegrated, label = FALSE, repel = TRUE, reduction = "umap",group.by = categary) + theme_bw() +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x='UMAP 1', y='UMAP 2') +
 scale_color_manual(values=categary.color.pallet)
ggsave(file=paste0(categary,'.UMAP.jpg'),width=5.6,height=4.8,unit='in',dpi = 800)

Idents(STobj.GEXintegrated) <- 'GEXintegrated_snn_res.0.4'
DefaultAssay(STobj.GEXintegrated) <- 'SCT'
#Identification of Spatially Variable Features
#Seurat offers two workflows to identify molecular features that correlate with spatial location within a tissue.
#The first is to perform differential expression based on pre-annotated anatomical regions within the tissue,
#which may be determined either from unsupervised clustering or prior knowledge.
#This strategy works will in this case, as the clusters above exhibit clear spatial restriction.
de_markers <- FindMarkers(STobj.GEXintegrated, ident.1 = 1, ident.2 = 2,test.use = "wilcox")
c0_markers <- FindMarkers(object = STobj.GEXintegrated, ident.1 = 0,only.pos =T)
c1_markers <- FindMarkers(object = STobj.GEXintegrated, ident.1 = 1,only.pos =T)
c2_markers <- FindMarkers(object = STobj.GEXintegrated, ident.1 = 2,only.pos =T)
c3_markers <- FindMarkers(object = STobj.GEXintegrated, ident.1 = 3,only.pos =T)
c4_markers <- FindMarkers(object = STobj.GEXintegrated, ident.1 = 4,only.pos =T)
c5_markers <- FindMarkers(object = STobj.GEXintegrated, ident.1 = 5,only.pos =T)
c6_markers <- FindMarkers(object = STobj.GEXintegrated, ident.1 = 6,only.pos =T)
c7_markers <- FindMarkers(object = STobj.GEXintegrated, ident.1 = 7,only.pos =T)
c8_markers <- FindMarkers(object = STobj.GEXintegrated, ident.1 = 8,only.pos =T)


c0_markers$marker <- rownames(c0_markers)
c1_markers$marker <- rownames(c1_markers)
c2_markers$marker <- rownames(c2_markers)
c3_markers$marker <- rownames(c3_markers)
c4_markers$marker <- rownames(c4_markers)
c5_markers$marker <- rownames(c5_markers)
c6_markers$marker <- rownames(c6_markers)
c7_markers$marker <- rownames(c7_markers)
c8_markers$marker <- rownames(c8_markers)


c0_markers$cluster <- rep('c0',nrow(c0_markers))
c1_markers$cluster <- rep('c1',nrow(c1_markers))
c2_markers$cluster <- rep('c2',nrow(c2_markers))
c3_markers$cluster <- rep('c3',nrow(c3_markers))
c4_markers$cluster <- rep('c4',nrow(c4_markers))
c5_markers$cluster <- rep('c5',nrow(c5_markers))
c6_markers$cluster <- rep('c6',nrow(c6_markers))
c7_markers$cluster <- rep('c7',nrow(c7_markers))
c8_markers$cluster <- rep('c8',nrow(c8_markers))


allMarkers <- rbind(c0_markers,c1_markers,c2_markers,c3_markers,c4_markers,c5_markers,c6_markers,c7_markers,c8_markers)
write.table(allMarkers,file='allMarkers.tsv',sep='\t',quote=F,row.names=F,col.names=T)

#Gene expression visualization
DefaultAssay(STobj.GEXintegrated) <- 'SCT'
feature='CFH'
pdf(file=paste0(feature,'_spatialFearturePlot.alpha0.1.pdf'))
#pt.size.factor- This will scale the size of the spots. Default is 1.6
#alpha - minimum and maximum transparency. Default is c(1, 1).
SpatialFeaturePlot(STobj.GEXintegrated, features = feature, pt.size.factor = 1.6,alpha = c(0.1, 1) ) + theme(legend.position = "right")
dev.off()

cluster.color.pallet <- as.character(c('#1f77b4','#ff7f0e','#2ca02c','#d62728','#7733B7',
  '#8c564b','#E236AF','#8C8C00','#700B35','#06BACE',
  "#00D8FF", "#992756", "#DADD00", "#F76AA7", "#C4796A",
  "#246A73", "#FF6360", "#4DDD30", "#FFA449", "#74A6E8",
  "#955BC1","#FF0044","#F6757A","#265C42","#193C3E",
  "#63C74D"))

SpatialDimPlot(STobj.GEXintegrated, label = FALSE, cols=cluster.color.pallet)
ggsave(file="Spatial_cluster_DimPlot.pdf")

# Integration with single-cell data
load('/work/project/xuanyu/project/MFSmultiome/Signac/Integration/iter5/step2/seuObj.step2.Rdata')
DefaultAssay(seuObj) <- 'GEXintegrated'
anchors <- FindTransferAnchors(reference = seuObj, query =STobj.GEXintegrated, normalization.method = "SCT",reference.assay ="GEXintegrated")
predictions.assay <- TransferData(anchorset = anchors, refdata = seuObj$subcluster, prediction.assay = TRUE,
    weight.reduction = STobj.GEXintegrated[["pca"]], dims = 1:30)
STobj.GEXintegrated[["predictions"]] <- predictions.assay
#Now we get prediction scores for each spot for each subpopulation.
DefaultAssay(STobj.GEXintegrated) <- "predictions"
> rownames( STobj.GEXintegrated[["predictions"]])
 [1] "VSMC1"      "VSMC2"      "venousEC"   "FB1"        "lEC"
 [6] "VSMC3"      "Mac1"       "Mural1"     "Mural2"     "Tcell"
[11] "arteryEC"   "Mac2"       "capllaryEC" "FB2"        "Adipocyte"
[16] "FB3"        "NK"         "max"

feature <- 'NK'
SpatialPlot(object = STobj.GEXintegrated, features = feature, ncol = 2,alpha = c(0.1, 1))
ggsave(file=paste0(feature,'.spatialPlot.pdf'))

save(STobj.GEXintegrated,file='STobj.GEXintegrated.Rdata')
