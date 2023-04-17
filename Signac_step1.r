library(Signac)
library(Seurat)
library(future)
plan("multicore", workers = 10)
options(future.globals.maxSize = 200 * 1024 ^ 3)
plan()

load('/work/project/xuanyu/project/MFSmultiome/Signac/creatObj/seuObj.combined.original.Rdata')

load('/work/project/xuanyu/project/MFSmultiome/Signac/QC/AD1565/AD1565_cells.Rdata')
load('/work/project/xuanyu/project/MFSmultiome/Signac/QC/AD1961/AD1961_cells.Rdata')
load('/work/project/xuanyu/project/MFSmultiome/Signac/QC/AD2069/AD2069_cells.Rdata')
load('/work/project/xuanyu/project/MFSmultiome/Signac/QC/AD2099/AD2099_cells.Rdata')
load('/work/project/xuanyu/project/MFSmultiome/Signac/QC/YZ18/YZ18_cells.Rdata')
load('/work/project/xuanyu/project/MFSmultiome/Signac/QC/YZ19/YZ19_cells.Rdata')
load('/work/project/xuanyu/project/MFSmultiome/Signac/QC/YZ30/YZ30_cells.Rdata')
load('/work/project/xuanyu/project/MFSmultiome/Signac/QC/YZ7/YZ7_cells.Rdata')

allCleanCells <- c(AD1565_cells,AD1961_cells,AD2069_cells,AD2099_cells,YZ18_cells,YZ19_cells,YZ30_cells,YZ7_cells)

VSMC_cells <- get(load('/work/project/xuanyu/project/MFSmultiome/Signac/Integration/subclustering/VSMC/latest/VSMC_cleanCell.Rdata'))
FB_cells <- get(load('/work/project/xuanyu/project/MFSmultiome/Signac/Integration/subclustering/FB/FBcells.clean.Rdata'))
Myeloid_cells <- get(load('/work/project/xuanyu/project/MFSmultiome/Signac/Integration/subclustering/Myeloid/MyleoidCells.Rdata'))
Lymphoid_cells <- get(load('/work/project/xuanyu/project/MFSmultiome/Signac/Integration/subclustering/Lymphoid/LymphoidCells.Rdata'))
vEC_cells <- get(load('/work/project/xuanyu/project/MFSmultiome/Signac/Integration/subclustering/vEC/vECCells.Rdata'))
lEC_cells <- get(load('/work/project/xuanyu/project/MFSmultiome/Signac/Integration/subclustering/lEC/lECCells.Rdata'))
Adipocyte_cells <- get(load('/work/project/xuanyu/project/MFSmultiome/Signac/Integration/subclustering/Adipocyte/AdipocyteCells.Rdata'))
Mural_cells <- get(load('/work/project/xuanyu/project/MFSmultiome/Signac/Integration/subclustering/mural/muralCells.Rdata'))

cleanCells <- c(VSMC_cells,FB_cells,Myeloid_cells,Lymphoid_cells,vEC_cells,lEC_cells,Adipocyte_cells,Mural_cells)

allCleanCells <- allCleanCells[allCleanCells %in% cleanCells]

seuObj.combined@meta.data$clean <- seuObj.combined@meta.data$CB %in% allCleanCells
seuObj.combined <- subset(seuObj.combined, clean=="TRUE")
# split the dataset into a list of seurat objects
seuObj.list <- SplitObject(seuObj.combined, split.by = "orig.ident")

# ---------- Gene expression data processing ------------------
DefaultAssay(seuObj.combined) <- "RNA"
# Normalize datasets individually by SCTransform(), instead of NormalizeData() prior to integration
# SCTransform (Hafemeister and Satija, 2019) procedure omits the need for heuristic steps including pseudocount addition or log-transformation
# The latest version of sctransform also supports using glmGamPoi package which substantially improves the speed of the learning procedure
#ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
# SCTransform is designed to eliminate the influence of sequencing depth. The nCount_RNA is used to construct the sct mode, so you don't need to put it into vars.to.regress.
seuObj.list <- lapply(X = seuObj.list, function(i){SCTransform(i,vars.to.regress = c("percent.mito","S.Score", "G2M.Score"), method = "glmGamPoi",assay='RNA')})
# Typically use 3,000 or more features for analysis downstream of sctransform
features <- SelectIntegrationFeatures(object.list = seuObj.list, nfeatures = 3000, assay =c('SCT','SCT','SCT','SCT','SCT','SCT','SCT','SCT'))
# Run the PrepSCTIntegration() function prior to identifying anchors
seuObj.list <- PrepSCTIntegration(object.list = seuObj.list, anchor.features = features, assay =c('SCT','SCT','SCT','SCT','SCT','SCT','SCT','SCT'))

# integration
GEX.anchors <- FindIntegrationAnchors(object.list = seuObj.list, normalization.method = "SCT", anchor.features = features)
seuObj.GEXintegrated <- IntegrateData(anchorset = GEX.anchors, normalization.method = "SCT",new.assay.name = "GEXintegrated",dims = 1:30)

#Perform dimensionality reduction by PCA and UMAP embedding
seuObj.GEXintegrated <- RunPCA(seuObj.GEXintegrated, verbose = FALSE, assay = "GEXintegrated",npcs = 50)


# ----------------- DNA accessibility data processing ---------------------
# compute LSI
# 'q5' to set the top 95% most common features as the VariableFeatures; setting to 10 will include features in >10 cells
# frequency-inverse document frequency (TF-IDF) normalization. This is a two-step normalization procedure,
# that both normalizes across cells to correct for differences in cellular sequencing depth,
# and across peaks to give higher values to more rare peaks.
# run singular value decomposition (SVD) on the TD-IDF matrix,
# Feature selection; Normalization, Dimension reduction
# The combined steps of TF-IDF followed by SVD are known as latent semantic indexing (LSI)

# # process the combined dataset
DefaultAssay(seuObj.combined) <- "ATAC"

seuObj.combined <- FindTopFeatures(seuObj.combined, min.cutoff = 'q5')
seuObj.combined <- RunTFIDF(seuObj.combined)
seuObj.combined <- RunSVD(seuObj.combined)
#seuObj.combined <- RunUMAP(seuObj.combined, reduction = "lsi", dims = 2:30)
#DimPlot(seuObj.combined, label = TRUE, repel = TRUE, reduction = "umap",group.by = "orig.ident")

DefaultAssay(seuObj.list$AD1565) <- "ATAC"
seuObj.list$AD1565 <- FindTopFeatures(seuObj.list$AD1565, min.cutoff = 'q5'); seuObj.list$AD1565 <- RunTFIDF(seuObj.list$AD1565); seuObj.list$AD1565 <- RunSVD(seuObj.list$AD1565)
DefaultAssay(seuObj.list$AD1961) <- "ATAC"
seuObj.list$AD1961 <- FindTopFeatures(seuObj.list$AD1961, min.cutoff = 'q5'); seuObj.list$AD1961 <- RunTFIDF(seuObj.list$AD1961); seuObj.list$AD1961 <- RunSVD(seuObj.list$AD1961)
DefaultAssay(seuObj.list$AD2069) <- "ATAC"
seuObj.list$AD2069 <- FindTopFeatures(seuObj.list$AD2069, min.cutoff = 'q5'); seuObj.list$AD2069 <- RunTFIDF(seuObj.list$AD2069); seuObj.list$AD2069 <- RunSVD(seuObj.list$AD2069)
DefaultAssay(seuObj.list$AD2099) <- "ATAC"
seuObj.list$AD2099 <- FindTopFeatures(seuObj.list$AD2099, min.cutoff = 'q5'); seuObj.list$AD2099 <- RunTFIDF(seuObj.list$AD2099); seuObj.list$AD2099 <- RunSVD(seuObj.list$AD2099)
DefaultAssay(seuObj.list$YZ18) <- "ATAC"
seuObj.list$YZ18 <- FindTopFeatures(seuObj.list$YZ18, min.cutoff = 'q5'); seuObj.list$YZ18 <- RunTFIDF(seuObj.list$YZ18); seuObj.list$YZ18 <- RunSVD(seuObj.list$YZ18)
DefaultAssay(seuObj.list$YZ19) <- "ATAC"
seuObj.list$YZ19 <- FindTopFeatures(seuObj.list$YZ19, min.cutoff = 'q5'); seuObj.list$YZ19 <- RunTFIDF(seuObj.list$YZ19); seuObj.list$YZ19 <- RunSVD(seuObj.list$YZ19)
DefaultAssay(seuObj.list$YZ30) <- "ATAC"
seuObj.list$YZ30 <- FindTopFeatures(seuObj.list$YZ30, min.cutoff = 'q5'); seuObj.list$YZ30 <- RunTFIDF(seuObj.list$YZ30); seuObj.list$YZ30 <- RunSVD(seuObj.list$YZ30)
DefaultAssay(seuObj.list$YZ7) <- "ATAC"
seuObj.list$YZ7 <- FindTopFeatures(seuObj.list$YZ7, min.cutoff = 'q5'); seuObj.list$YZ7 <- RunTFIDF(seuObj.list$YZ7); seuObj.list$YZ7 <- RunSVD(seuObj.list$YZ7)

# find integration anchors using reciprocal LSI projection
# An important first step in any integrative analysis of single-cell chromatin data is to ensure that the same features are measured in each dataset.

integration.anchors <- FindIntegrationAnchors(
  object.list = seuObj.list,
  anchor.features = rownames(seuObj.list[[1]]),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
seuObj.ATACintegrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = seuObj.combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)


# Compute the correlation between total counts and each reduced dimension component.
# The first LSI component often captures sequencing depth (technical variation) rather than biological variation. If this is the case, the component should be removed from downstream analysis
DepthCor(seuObj.ATACintegrated, assay = 'ATAC', reduction = "integrated_lsi", n = 50)
# create a new UMAP using the integrated embeddings
#seuObj.ATACintegrated.UMAP <- RunUMAP(seuObj.ATACintegrated, reduction = "integrated_lsi", dims = 2:30)


# Joint UMAP visualization
# Using the weighted nearest neighbor methods, we can compute a joint neighbor graph that represent
# both the gene expression and DNA accessibility measurements.
# build a joint neighbor graph using both assays

seuObj <- seuObj.GEXintegrated
seuObj@assays$ATAC <- seuObj.ATACintegrated@assays$ATAC
seuObj@reductions$lsi <- seuObj.combined@reductions$lsi
seuObj@reductions$integrated_lsi <- seuObj.ATACintegrated@reductions$integrated_lsi

seuObj <- FindMultiModalNeighbors(
  object = seuObj,
  reduction.list = list("pca", "integrated_lsi"),
  dims.list = list(1:30, 2:30),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

DefaultAssay(seuObj) <- 'SCT'
# build a joint UMAP visualization
seuObj <- RunUMAP(
  object = seuObj,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

# clustering based on GEX data
seuObj <- FindNeighbors(object = seuObj, reduction = 'pca', dims = 1:30)
seuObj <- FindClusters(object = seuObj, verbose = FALSE, algorithm = 1,graph.name='GEXintegrated_snn',resolution =1)
seuObj <- FindClusters(object = seuObj, verbose = FALSE, algorithm = 1,graph.name='GEXintegrated_snn',resolution =0.6)
seuObj <- FindClusters(object = seuObj, verbose = FALSE, algorithm = 1,graph.name='GEXintegrated_snn',resolution =0.4)
seuObj <- FindClusters(object = seuObj, verbose = FALSE, algorithm = 1,graph.name='GEXintegrated_snn',resolution =0.2)
seuObj <- FindClusters(object = seuObj, verbose = FALSE, algorithm = 1,graph.name='GEXintegrated_snn',resolution =0.1)

#map celltype
library(plyr)
cleanCells <- c(VSMC_cells,FB_cells,Myeloid_cells,Lymphoid_cells,vEC_cells,lEC_cells,Adipocyte_cells,Mural_cells)
seuObj@meta.data$cellType <- mapvalues(seuObj@meta.data$CB,from=cleanCells,to=c(rep('VSMC',length(VSMC_cells)),rep('FB',length(FB_cells)),
rep('Myeloid',length(Myeloid_cells)),rep('Lymphoid',length(Lymphoid_cells)),rep('vEC',length(vEC_cells)),
rep('lEC',length(lEC_cells)),rep('Adipocyte',length(Adipocyte_cells)),rep('Mural',length(Mural_cells))))

load('/work/project/xuanyu/project/MFSmultiome/Signac/Integration/subclustering/VSMC/latest/VSMC_subcluster_CB.Rdata')
load('/work/project/xuanyu/project/MFSmultiome/Signac/Integration/subclustering/FB/FB_subcluster_CB.Rdata')
load('/work/project/xuanyu/project/MFSmultiome/Signac/Integration/subclustering/Myeloid/Myeloid_subcluster_CB.Rdata')
load('/work/project/xuanyu/project/MFSmultiome/Signac/Integration/subclustering/Lymphoid/Lymphoid_subcluster_CB.Rdata')
load('/work/project/xuanyu/project/MFSmultiome/Signac/Integration/subclustering/mural/Mural_subcluster_CB.Rdata')
load('/work/project/xuanyu/project/MFSmultiome/Signac/Integration/subclustering/vEC/vEC_subcluster_CB.Rdata')
load('/work/project/xuanyu/project/MFSmultiome/Signac/Integration/subclustering/Adipocyte/Adipocyte_subcluster_CB.Rdata')
load('/work/project/xuanyu/project/MFSmultiome/Signac/Integration/subclustering/lEC/lEC_subcluster_CB.Rdata')

subcluster <- c(VSMC_subcluster,FB_subcluster,Myeloid_subcluster,Lymphoid_subcluster,Mural_subcluster,vEC_subcluster,Adipocyte_subcluster,lEC_subcluster)
CB <- c(VSMC_CB,FB_CB,Myeloid_CB,Lymphoid_CB,Mural_CB,vEC_CB,Adipocyte_CB,lEC_CB)
seuObj@meta.data$subcluster <- mapvalues(seuObj@meta.data$CB,from=CB,to=subcluster)

sub <-c("Adipocyte","arteryEC","capllaryEC","FB1","FB2",
"FB3","lEC","Mac1","Mac2","Mural1",
"Mural2","NK","Tcell","venousEC","VSMC1",
"VSMC2","VSMC3")
cellType2 <- c("Adipocyte","vEC","vEC","FB","FB",
"FB","lEC","Macrophage","Macrophage","Mural",
"Mural","NK","T cell","vEC","VSMC",
"VSMC","VSMC")
seuObj@meta.data$cellType2 <- mapvalues(seuObj@meta.data$subcluster,from=sub,to=cellType2)


# visualization
# Dimplot
library(ggplot2)

# UMAP categary data plot
categary='subcluster'
categary.color.pallet  <- as.character(c('#1f77b4','#ff7f0e','#2ca02c','#d62728','#7733B7',
  '#8c564b','#E236AF','#8C8C00','#700B35','#06BACE',
"#193C3E", "#992756", "#DADD00", "#F76AA7", "#C4796A",
  "#246A73", "#FF6360", "#4DDD30", "#FFA449", "#74A6E8",
  "#955BC1","#FF0044","#F6757A","#265C42","#00D8FF",
  "#63C74D","#800026","#bd0026"))
DimPlot(seuObj, label = FALSE, repel = TRUE, reduction = "umap",group.by = categary) + theme_bw() +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x='UMAP 1', y='UMAP 2') +
 scale_color_manual(values=categary.color.pallet)
ggsave(file='subcluster.UMAP.jpg',width=5.6,height=4.8,unit='in',dpi = 600)
categary='cellType'
DimPlot(seuObj, label = TRUE, repel = TRUE, reduction = "umap",group.by = categary) + theme_bw() +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x='UMAP 1', y='UMAP 2') +
 scale_color_manual(values=categary.color.pallet)
ggsave(file='cellType.UMAP.jpg',width=5.6,height=4.8,unit='in',dpi = 600)
categary='cellType2'
DimPlot(seuObj, label = FALSE, repel = TRUE, reduction = "umap",group.by = categary) + theme_bw() +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x='UMAP 1', y='UMAP 2') +
 scale_color_manual(values=categary.color.pallet)
ggsave(file='cellType2.UMAP.jpg',width=5.6,height=4.8,unit='in',dpi = 600)

# set default assay and idents
DefaultAssay(seuObj) <- "SCT"
seuObj@meta.data$cellType <- factor(seuObj@meta.data$cellType,
levels=c('VSMC','FB','vEC','Myeloid','Lymphoid','Mural','Adipocyte','lEC'))
Idents(seuObj) <- "cellType"


# feature plot
feature='RYR2'
feature <- 'percent.mito'
feature <- 'nCount_ATAC'
feature <- 'nCount_RNA'
feature <- 'TSS.enrichment'
FeaturePlot(seuObj, features = feature, reduction = 'umap',cols =c('grey','red'))  +
 ggtitle(feature) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x='UMAP 1', y='UMAP 2')
ggsave(file=paste0(feature,'.UMAP.jpg'),width=5.6,height=4.82,unit='in',dpi = 600)

# Visualize co-expression of two features simultaneously
FeaturePlot(seuObj, features = c("DCN", "MYH11"), blend = TRUE)

# find all markers
markers.all <- FindAllMarkers(object = seuObj, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5, test.use="bimod")
markers <- FindMarkers(object = seuObj, ident.1='diseased',ident.2='healthy',only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use="bimod")
head(subset(markers.all,cluster == '14'),n=20)
save(markers.all,file='markers.all.cellType2.Rdata')
write.table(markers.all,file="markers.all.cellType2.tsv",sep="\t",row.names=F,quote=F,col.names=T)

# violin plot
#stack violin plot
features <- c('MYH11','DCN','VWF','ABCC9','ACACB','GPAM','ITGAM','NCR1','CD3D','MMRN1')
features <- c('MYH11','DCN','VWF','F13A1','SKAP1','ADGRL3','ACACB','MMRN1')
violin.color.pallet  <- as.character(c('#1f77b4','#ff7f0e','#2ca02c','#d62728','#7733B7',
  '#8c564b','#E236AF','#8C8C00','#700B35','#06BACE',
"#193C3E", "#992756", "#DADD00", "#F76AA7", "#C4796A",
  "#246A73", "#FF6360", "#4DDD30", "#FFA449", "#74A6E8",
  "#955BC1","#FF0044","#F6757A","#265C42","#00D8FF",
  "#63C74D","#800026","#bd0026"))
VlnPlot(seuObj, assay = 'SCT', features = features,stack= TRUE,flip=TRUE) + theme_bw() +
theme(legend.position='none') +  scale_fill_manual(values=violin.color.pallet) + labs(x='',y='Normalized Gene Expression')
ggsave(file='cellType_marker_geneExp_violinPlot.pdf')
VlnPlot(seuObj, features = feature, split.by = "groups")


# dot plot
DotPlot(seuObj, features = features,cols =c('grey','red')) + RotatedAxis() + ggmin::theme_powerpoint() +labs(x='',y='')+
 theme(axis.text.x = element_text(colour = "black",angle = 30, size = 10,vjust=0.5))
ggsave(file='cellType_marker_geneExp_DotPlot.pdf')
# SplitDotPlotGG has been replaced with the `split.by` parameter for DotPlot
DotPlot(seuObj, features = features) + RotatedAxis() + ggmin::theme_powerpoint()
DotPlot(seuObj, features = features, split.by = "group") + RotatedAxis() + ggmin::theme_powerpoint()

# Interactive plotting features
plot <- FeaturePlot(pbmc3k.final, features = "MS4A1")
HoverLocator(plot = plot, information = FetchData(pbmc3k.final, vars = c("ident", "PC_1", "nFeature_RNA")))


# save object
save(seuObj,file='seuObj.step1.Rdata')

# UMAP embedding soly based on GEX data
seuObj.GEXintegrated <- RunUMAP(seuObj.GEXintegrated, reduction = "pca", assay = "GEXintegrated", dims = 1:30)
seuObj.GEXintegrated <- FindNeighbors(seuObj.GEXintegrated, reduction = "pca", dims = 1:30)
seuObj.GEXintegrated <- FindClusters(seuObj.GEXintegrated, resolution = 1)
seuObj.GEXintegrated@meta.data$cellType <- seuObj@meta.data$cellType
seuObj.GEXintegrated@meta.data$subcluster <- seuObj@meta.data$subcluster
seuObj.GEXintegrated@meta.data$cellType2 <- seuObj@meta.data$cellType2
DimPlot(seuObj.GEXintegrated, label = TRUE, repel = TRUE, reduction = "umap",group.by = "subcluster") + theme_bw() +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x='UMAP 1', y='UMAP 2')+
 scale_color_manual(values=categary.color.pallet)
ggsave(file='GEX_subcluster.UMAP.jpg',width=5.6,height=4.8,unit='in',dpi = 600)
DimPlot(seuObj.GEXintegrated, label = TRUE, repel = TRUE, reduction = "umap",group.by = "cellType") + theme_bw() +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x='UMAP 1', y='UMAP 2')+
 scale_color_manual(values=categary.color.pallet)
ggsave(file='GEX_cellType.UMAP.jpg',width=5.6,height=4.8,unit='in',dpi = 600)
# save
save(seuObj.GEXintegrated,file='seuObj.GEXintegrated.Rdata')

# UMAP embedding soly based on ATAC data
seuObj.ATACintegrated <- RunUMAP(seuObj.ATACintegrated, reduction = "integrated_lsi", dims = 2:30)
seuObj.ATACintegrated <- FindNeighbors(object = seuObj.ATACintegrated, reduction = 'integrated_lsi', dims = 2:30)
seuObj.ATACintegrated <- FindClusters(object = seuObj.ATACintegrated, verbose = FALSE, algorithm = 1,graph.name='ATAC_snn',
resolution =1)
seuObj.ATACintegrated@meta.data$cellType <- seuObj@meta.data$cellType
seuObj.ATACintegrated@meta.data$cellType2 <- seuObj@meta.data$cellType2
seuObj.ATACintegrated@meta.data$subcluster <- seuObj@meta.data$subcluster
DimPlot(seuObj.ATACintegrated, label = TRUE, repel = TRUE, reduction = "umap",group.by = "subcluster") + theme_bw() +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x='UMAP 1', y='UMAP 2')+
 scale_color_manual(values=categary.color.pallet)
ggsave(file='ATAC_subcluster.UMAP.jpg',width=5.6,height=4.8,unit='in',dpi = 600)
DimPlot(seuObj.ATACintegrated, label = TRUE, repel = TRUE, reduction = "umap",group.by = "cellType") + theme_bw() +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x='UMAP 1', y='UMAP 2')+
 scale_color_manual(values=categary.color.pallet)
ggsave(file='ATAC_cellType.UMAP.jpg',width=5.6,height=4.8,unit='in',dpi = 600)
# save
save(seuObj.ATACintegrated,file='seuObj.ATACintegrated.Rdata')
