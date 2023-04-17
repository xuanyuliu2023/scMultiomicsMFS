library(pheatmap)
library(AUCell)
df <- readRDS("int/3.4_regulonAUC.Rds")
load('/work/project/xuanyu/project/HCMscRNAseq/seurat/meta.data.table.RData')

CM1_cells <- row.names(subset(meta.data.table,subcluster=="CM1"))
CM2_cells <- row.names(subset(meta.data.table,subcluster=="CM2"))

allCellId <- c(CM1_cells,CM2_cells)
allCellId <- allCellId[allCellId %in% colnames(df)]
df <- df[,allCellId]

library(plyr)
clusterInfo <- mapvalues(allCellId, from=meta.data.table$cellID, to=meta.data.table$subcluster)
annoCol <- data.frame(cluster=clusterInfo)
row.names(annoCol)<- allCellId

library(RColorBrewer)
cc = colorRampPalette(rev(brewer.pal(n = 7,
     name = "RdYlBu")))



cluster.parameters=list(aucScore_matrix,scale = "row",cellwidth = 0.02, color = colorRampPalette(colors = c("blue","white","red"))(50),
                        cellheight = 6,fontsize_row=8,cluster_cols = F,fontsize_col=8,show_colnames=F,
                        cluster_rows = T,annotation_col=annoCol,display_numbers = F)
do.call("pheatmap", cluster.parameters)
do.call("pheatmap", c(cluster.parameters, filename='DEG_heatmap.pdf'))

-----------
# heatmap for Average Regulon Activity
library(SCENIC)
library(AUCell)
library(pheatmap)
scenicOptions <- readRDS(file="int/scenicOptions.Rds")
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
load('/work/project/xuanyu/project/HCMscRNAseq/seurat/meta.data.table.RData')
FB1_cells <- row.names(subset(meta.data.table,subcluster=="FB1"))
FB2_cells <- row.names(subset(meta.data.table,subcluster=="FB2"))
FB3_cells <- row.names(subset(meta.data.table,subcluster=="FB3"))
FB4_cells <- row.names(subset(meta.data.table,subcluster=="FB4"))
allCellId <- c(FB1_cells,FB2_cells,FB3_cells,FB4_cells)
allCellId <- allCellId[allCellId %in% colnames(regulonAUC)]
library(plyr)
clusterInfo <- mapvalues(allCellId, from=meta.data.table$cellID, to=meta.data.table$subcluster)
groupInfo <- mapvalues(allCellId, from=meta.data.table$cellID, to=meta.data.table$group)
annoCol <- data.frame(cluster=clusterInfo,group=groupInfo)
row.names(annoCol)<- allCellId

regulonActivity_byCellType <- sapply(split(rownames(annoCol), annoCol$cluster),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))

TFs <- row.names(regulonActivity_byCellType)
highCofReglons <- TFs[!grepl("_extended",TFs)]
regulonActivity_byCellType_highConf <- regulonActivity_byCellType[highCofReglons ,]

regulonActivity_byCellType_highConf_Scaled <- t(scale(t(regulonActivity_byCellType_highConf), center = T, scale=T))
regulonActivity_byCellType_highConf_Scaled_logic <- regulonActivity_byCellType_highConf_Scaled >0
regulonActivity_byCellType_highConf_Scaled_logic_specifc <- apply(regulonActivity_byCellType_highConf_Scaled_logic,1,sum)==1
regulonActivity_byCellType_highConf_Scaled.specifc <- regulonActivity_byCellType_highConf_Scaled[regulonActivity_byCellType_highConf_Scaled_logic_specifc,]
#regulonActivity_byCellType_Scaled <- regulonActivity_byCellType_Scaled[,c("c3","c1","c0","c5")]
#,breaks=seq(-0.01, 0.01, length.out = 100),
cluster.parameters <- list(regulonActivity_byCellType_highConf_Scaled.specific, fontsize_row=6,
                   color=colorRampPalette(c("blue","white","red"))(100),
                   treeheight_row=20, treeheight_col=20, border_color=NA, cellwidth = 20, cellheight = 10, cluster_cols = F, cluster_rows = T)
do.call("pheatmap", cluster.parameters)

cluster.parameters <- list(regulonActivity_byCellType_highConf_Scaled, fontsize_row=6,
                   color=colorRampPalette(c("blue","white","red"))(100),
                   treeheight_row=20, treeheight_col=20, border_color=NA, cellwidth = 20, cellheight = 6, cluster_cols = F, cluster_rows = T)
do.call("pheatmap", cluster.parameters)


addInfo <- regulonActivity_byCellType_highConf_Scaled[c("GABPA (14g)","ETV5 (47g)","HES1 (10g)","HLF (36g)","GATA4 (11g)","TCF7L2 (21g)"),]
regulonActivity_byCellType_highConf_Scaled.specific.final <- rbind(regulonActivity_byCellType_highConf_Scaled.specifc,addInfo)

regulons <- c("GABPA (14g)","ETV5 (47g)","HES1 (10g)","HLF (36g)","TCF7L2 (21g)","GATA4 (11g)",
"ZNF91 (16g)", "CREB3L1 (277g)","USF2 (17g)","CREB3L2 (47g)","LEF1 (143g)","MAX (42g)","PHF8 (245g)","ZNF454 (49g)","NRF1 (78g)","PRRX2 (98g)","NFKB2 (87g)",
"ZNF135 (13g)","POLR2A (1404g)","CEBPG (40g)","TAF1 (258g)","IRF2 (126g)","IRF1 (82g)", "NR3C1 (1927g)",
"JUNB (23g)","TEAD1 (10g)", "STAT3 (55g)", "ZNF75D (43g)", "HOMEZ (11g)", "HLTF (15g)", "KLF2 (119g)", "JUND (13g)", "SREBF1 (26g)","PAX8 (10g)")
regulonActivity_byCellType_highConf_Scaled.specific.final <- regulonActivity_byCellType_highConf_Scaled.specific.final[regulons,]

cluster.parameters <- list(regulonActivity_byCellType_highConf_Scaled.specific.final, fontsize_row=6,
                   color=colorRampPalette(c("blue","white","red"))(100),
                   treeheight_row=20, treeheight_col=20, border_color=NA, cellwidth = 20, cellheight = 10, cluster_cols = F, cluster_rows = F)
do.call("pheatmap", cluster.parameters)
do.call("pheatmap", c(cluster.parameters, filename='averageRegulonActivityHeatmap.pdf'))
save(regulonActivity_byCellType_highConf_Scaled.specific.final,file="regulonActivity_byCellType_highConf_Scaled.specific.final.Rdata")



############ by group
regulonActivity_byGroup <- sapply(split(rownames(annoCol), annoCol$group),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))

TFs <- row.names(regulonActivity_byGroup)
highCofReglons <- TFs[!grepl("_extended",TFs)]
regulonActivity_byGroup_highConf <- regulonActivity_byGroup[highCofReglons ,]
regulonActivity_byGroup_highConf <- as.data.frame(regulonActivity_byGroup_highConf)
regulonActivity_byGroup_highConf$diff <- regulonActivity_byGroup_highConf[,1]-regulonActivity_byGroup_highConf[,2]
regulonActivity_byGroup_highConf$fc <- log2(regulonActivity_byGroup_highConf[,1]/regulonActivity_byGroup_highConf[,2])
regulonActivity_byGroup_highConf <- regulonActivity_byGroup_highConf[order(regulonActivity_byGroup_highConf$diff,decreasing=T),]
regulonActivity_byGroup_highConf <- regulonActivity_byGroup_highConf[order(regulonActivity_byGroup_highConf$fc,decreasing=T),]
regulonActivity_byGroup_highConf.final <- subset(regulonActivity_byGroup_highConf,abs(fc)>0.5) #abs log2 fold change >0.5
regulonActivity_byGroup_highConf.final.genes <- unlist(strsplit(row.names(regulonActivity_byGroup_highConf.final)," "))[seq(1,82,2)]

DEGs <- read.table(file='/work/project/xuanyu/project/HCMscRNAseq/subClustering/FB/diffExpr/results.classified.simplified.tsv',sep='\t',header=T)
DEGinfo <- subset(DEGs, gene %in% regulonActivity_byGroup_highConf.final.genes)
row.names(DEGinfo) <- DEGinfo$gene
DEGinfo <- DEGinfo[regulonActivity_byGroup_highConf.final.genes,]
combindInfo <- cbind(regulonActivity_byGroup_highConf.final,DEGinfo)
combindInfo <- rbind(subset(combindInfo,fc > 0 & log2FC >1 & State == "up"),subset(combindInfo,fc < 0 & log2FC < -1 & State == "down"))
combindInfo$fc <- 0-combindInfo$fc
combindInfo.concise <- combindInfo[,c('fc','log2FC')]
combindInfo.concise$reglon <- row.names(combindInfo.concise)

library(ggplot2)
library(reshape2)
data.plot <- melt(combindInfo.concise)
data.plot$reglon <- factor(data.plot$reglon,levels=rev(combindInfo.concise$reglon))
ggplot(data.plot,aes(x=reglon,y=value,fill=variable))+geom_bar(stat="identity") +theme_bw() + scale_y_continuous(breaks = seq(-2, 6, 0.5)) + coord_flip()
ggsave(file="reglonActivityDEG.barplot.pdf")
