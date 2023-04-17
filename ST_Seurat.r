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

seuObj.list$AD2242@meta.data$keep <-  paste0(colnames(seuObj.list$AD2242),'_1') %in% spots_SC0_SC1_SC4
seuObj.list$AD2242 <- subset(seuObj.list$AD2242, keep == TRUE)

seuObj.list$AD2287@meta.data$keep <-  paste0(colnames(seuObj.list$AD2287),'_2') %in% spots_SC0_SC1_SC4
seuObj.list$AD2287 <- subset(seuObj.list$AD2287, keep == TRUE)

seuObj.list$YZ75@meta.data$keep <-  paste0(colnames(seuObj.list$YZ75),'_3') %in% spots_SC0_SC1_SC4
seuObj.list$YZ75 <- subset(seuObj.list$YZ75, keep == TRUE)

seuObj.list$YZ74@meta.data$keep <-  paste0(colnames(seuObj.list$YZ74),'_4') %in% spots_SC0_SC1_SC4
seuObj.list$YZ74 <- subset(seuObj.list$YZ74, keep == TRUE)

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

library(plyr)
STobj.GEXintegrated@meta.data$seurat_clusters <- mapvalues(rownames(STobj.GEXintegrated@meta.data),from=rownames(meta.data.NormDat),to=meta.data.NormDat$seurat_clusters)
# UMAP categary data plot

categary.color.pallet  <- as.character(c('#1f77b4','#ff7f0e','#7733B7',
  '#8c564b','#E236AF','#8C8C00','#700B35','#06BACE',
"#193C3E", "#992756", "#DADD00", "#F76AA7", "#C4796A",
  "#246A73", "#FF6360", "#4DDD30", "#FFA449", "#74A6E8",
  "#955BC1","#FF0044","#F6757A","#265C42","#00D8FF",
  "#63C74D","#800026","#bd0026"))
categary='GEXintegrated_snn_res.0.2'
categary='sampleID'
categary='seurat_clusters'
DimPlot(STobj.GEXintegrated, label = FALSE, repel = TRUE, reduction = "umap",group.by = categary,pt.size=1) + theme_bw() +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x='UMAP 1', y='UMAP 2') +
 scale_color_manual(values=categary.color.pallet)
ggsave(file=paste0(categary,'.UMAP.jpg'),width=5.6,height=4.8,unit='in',dpi = 800)
save(STobj.GEXintegrated,file='STobj.GEXintegrated_SC0SC1SC4.RData')

---------------------------------------------
conda activate Monocle3
library(Seurat)
library(monocle3)
STdata <- get(load('STobj.GEXintegrated_SC0SC1SC4.RData'))

library(SingleCellExperiment)
DefaultAssay(STdata) <- 'SCT'

cell_metadata <- STdata@meta.data
expression_matrix<-STdata@assays$SCT@counts

gene_annotation <- data.frame(gene_short_name=rownames(expression_matrix))
rownames(gene_annotation) <- gene_annotation$gene_short_name
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

colData(cds)
rowData(cds)

# only consider a subset of clusters for identifying genes that change as a function of paseudotime
#cds_subset <- cds[,colData(cds)$seurat_clusters %in% c('SC0','SC1','SC4')]
cds_subset <- cds
# Pre-process the data
cds_subset <- preprocess_cds(cds_subset, num_dim = 30)
# Reduce dimensionality and visualize the cells
cds_subset <- reduce_dimension(cds_subset, reduction_method = 'UMAP')
# 指定UMAP
cds_subset@int_colData$reducedDims@listData$UMAP <- STdata@reductions$umap@cell.embeddings

#cds_subset@int_colData@listData$reducedDims$UMAP <- STdata@reductions$umap@cell.embeddings
plot_cells(cds_subset)
library(ggplot2)
plot_cells(cds_subset, color_cells_by="seurat_clusters") + scale_color_manual(values=c('#1f77b4','#ff7f0e','#7733B7','#8c564b','#E236AF','#8C8C00','#700B35')) + theme_bw()
#Group cells into clusters
#we run cluster_cells(), each cell is assigned not only to a cluster but also to a partition.
#When you are learning trajectories, each partition will eventually become a separate trajectory.
cds_subset = cluster_cells(cds_subset, reduction_method = 'UMAP', resolution=NULL)  # Note: the trajectory is dependent on cluster resolution
plot_cells(cds_subset, color_cells_by = "partition")
#Learn the trajectory graph
#To reduce the "branchyness" of your trajectories, you can adjust the ncenter parameter
cds_subset <- learn_graph(cds_subset,learn_graph_control=list(ncenter=100))

library(ggplot2)
pdf(file='seurat_clusters_trajectory_graph.pdf')
plot_cells(cds_subset,
           color_cells_by = "seurat_clusters",
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE,
           cell_size = 1,
           label_cell_groups = FALSE,
           trajectory_graph_segment_size = 1.5,
           label_branch_points=TRUE) + scale_color_manual(values=c('#1f77b4','#ff7f0e','#7733B7','#8c564b','#E236AF','#8C8C00','#700B35')) + theme_bw()
dev.off()

# In order to place the cells in order,
#we need to tell Monocle where the "beginning" of the biological process is.
# We do so by choosing regions of the graph that we mark as "roots" of the trajectory
#  In general, you should choose at least one root per partition.
# Order the cells in pseudotime
cds_subset <- order_cells(cds_subset)

pdf(file='pseudotime_trajectory_graph.pdf')
plot_cells(cds_subset,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           cell_size = 1,
           trajectory_graph_segment_size = 1.5) + theme_bw()
dev.off()

plot_cells(cds_subset,
           color_cells_by = "RYR2",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           cell_size = 1,
           trajectory_graph_segment_size = 1.5) + theme_bw()

# ---------------detect genes change along the pseudotime
cds_subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)
dim(cds_subset_pr_test_res)
library(dplyr)
cds_subset_pr_test_res <- cds_subset_pr_test_res %>% arrange(q_value,desc(morans_I))

# # only consider protein coding genes
geneSymbol2Type <- read.table(file="geneSymbol_geneType.tsv",sep='\t',header=F)
proteinGenes <- subset(geneSymbol2Type,V2=="protein_coding")$V1
cds_subset_pr_test_res  <- subset(cds_subset_pr_test_res, gene_short_name %in% proteinGenes)
sigGenes.protein.df <- subset(cds_subset_pr_test_res, q_value < 0.05) # threshold
#dim(sigGenes.protein.df)
write.table(cds_subset_pr_test_res,file='cds_subset_pr_test_res.tsv',quote=F,sep='\t',col.names = T,row.names = F)
write.table(sigGenes.protein.df,file='sigGenes.protein.df.tsv',quote=F,sep='\t',col.names = T,row.names = F)



# plot genes_trajectoryPlot
MARKER_genes <- c("RYR2","COL8A1","CDH11","TEAD1","BACH2","FOXN3","MYH11")
MARKER_genes <- c('CFH','CXCL12','AEBP1','NRXN3','TNFRSF11B','SERPINE1','BACH1')
MARKER_cds <- cds_subset[rowData(cds_subset)$gene_short_name %in% MARKER_genes,]
plot_genes_in_pseudotime(MARKER_cds,
                         color_cells_by="seurat_clusters",
                         cell_size = 0.1) + scale_color_manual(values=c('#1f77b4','#ff7f0e','#7733B7','#8c564b','#E236AF','#8C8C00','#700B35')) + theme_bw()
ggsave(file="plot_genes_in_pseudotime.pdf")
save(cds_subset,file="cds_subset.Rdata")
# save as  genes_trajectoryPlot.pdf size: 16 * 4
----------------------------
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3)
sigGenes.protein.df <- read.table(file='sigGenes.protein.df',header=T,sep='\t')


#modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
# <- row.names(subset(modulated_genes, q_value == 0 & morans_I > 0.25))
#genes
#genesTodiscard <- subset(out,Cluster=="cluster1")$GeneID
sigGenes.protein.df <- sigGenes.protein.df[sort(sigGenes.protein.df$morans_I,decreasing=T,index.return = T)$ix,]
#genes <- sigGenes.protein.df$gene_short_name[1:50] #only display the top 50 genes
genes <- sigGenes.protein.df$gene_short_name
#genes <- genes[!genes %in% genesTodiscard]

pt.matrix <- exprs(cds_subset)[match(genes,rownames(rowData(cds_subset))),order(pseudotime(cds_subset))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes;
#K means with 6 groups
set.seed(123)
htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = FALSE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 8),
  row_km = 6,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)
htkm <- draw(htkm)

pdf(file='htkm.heatmap.pdf')
print(htkm)
dev.off()

r.dend <- row_dend(htkm)  #Extract row dendrogram
rcl.list <- row_order(htkm)  #Extract clusters (output is a list)
# loop to extract genes for each cluster.
for (i in 1:length(row_order(htkm))){
  if (i == 1) {
    clu <- t(t(row.names(pt.matrix[row_order(htkm)[[i]],])))
    OUT <- cbind(clu, paste("cluster", i, sep=""))
    colnames(OUT) <- c("GeneID", "Cluster")
  } else {
    clu <- t(t(row.names(pt.matrix[row_order(htkm)[[i]],])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    OUT <- rbind(OUT, clu)
  }
}

out <- as.data.frame(OUT)
dim(out)
table(out$Cluster)
subset(out,Cluster=="cluster3")
subset(out,GeneID=="NPPB")
write.table(out,sep='\t',file='htkm_cluste_gene.tsv',quote=F,col.names = T,row.names = F)
--------------------

library(monocle3)
library(dplyr)
load('cds_subset.Rdata')
# the pseudotime values for all cells/barcode-spots are obtained via
pseudotime_vec <- monocle3::pseudotime(cds_subset)
pseudotime_df <-
  base::as.data.frame(pseudotime_vec)
# 2. create a barcodes-key variable and rename the pseudotime variable
feature_df <-
  magrittr::set_colnames(x = pseudotime_df, value = "Pseudotime") %>%
  tibble::rownames_to_column(var = "barcodes")

# 3. add to spata-object via
spataObj <- addFeatures(object = spataObj,
                        feature_names = "Pseudotime",
                        feature_df = feature_df,
                        overwrite = TRUE)
# visualize pseudotime on the surface
# note: the figure must be saved through Rstudio 'pseudotime_surfacePlot.pdf'
plotSurface(spataObj, color_by = "Pseudotime", smooth = TRUE, smooth_span = 0.2, pt_size = 2)

#----------
feature_df<- cbind(feature_df,colData(cds_subset))
feature_df_diseased <- subset(feature_df,group == "diseased")
ggplot(feature_df_diseased,aes(x=Pseudotime,fill=seurat_clusters)) + geom_density(alpha=0.5) + scale_fill_manual(values=c('#1f77b4','#ff7f0e','#7733B7','#8c564b','#E236AF','#8C8C00','#700B35')) + theme_bw() + labs(x='Pseudotime',y='Density')
feature_df_healthy <- subset(feature_df,group == "healthy")
ggplot(feature_df_healthy,aes(x=Pseudotime,fill=seurat_clusters)) + geom_density(alpha=0.5) + scale_fill_manual(values=c('#1f77b4','#ff7f0e','#7733B7','#8c564b','#E236AF','#8C8C00','#700B35')) + theme_bw() + labs(x='Pseudotime',y='Density')
ggplot(feature_df,aes(x=Pseudotime,fill=group)) + geom_density(alpha=0.5) + scale_fill_manual(values=c('#1f77b4','#ff7f0e','#7733B7','#8c564b','#E236AF','#8C8C00','#700B35')) + theme_bw() + labs(x='Pseudotime',y='Density')

ggplot(feature_df,aes(x=Pseudotime,fill=seurat_clusters)) + geom_density(alpha=0.5) + scale_fill_manual(values=c('#1f77b4','#ff7f0e','#7733B7','#8c564b','#E236AF','#8C8C00','#700B35')) + theme_bw() + labs(x='Pseudotime',y='Density')
ggsave(file='pseudotimeDensityPlot.cluster.pdf')
__________________________________________________
