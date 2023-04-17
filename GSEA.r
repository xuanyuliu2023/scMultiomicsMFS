library(Seurat)
library(plyr)
load('/work/project/xuanyu/project/MFSmultiome/Signac/Integration/iter5/step2/seuObj.step2.Rdata')

cellstobeconsidered.df <- subset(seuObj@meta.data,cellType == "VSMC")
cellstobeconsidered.df$cell <- row.names(cellstobeconsidered.df)

library(dplyr)
> table(cellstobeconsidered.df$group)

diseased  healthy
   13497    14449


#cellstobeconsidered.df.sampled.byGroup <- cellstobeconsidered.df %>% group_by(Group) %>% sample_n(1533)
#cellstobeconsidered <- cellstobeconsidered.df.sampled.byGroup$cell

HEALTHY.Group.cell <- subset(cellstobeconsidered.df,group=='healthy')$cell
DISEASED.Group.cell <- subset(cellstobeconsidered.df,group=='diseased')$cell
cellstobeconsidered  <- c(DISEASED.Group.cell,HEALTHY.Group.cell)
cellstobeconsidered.df <- cellstobeconsidered.df[cellstobeconsidered,]
sampleSize <- length(cellstobeconsidered)
sink("sample.phenotype.cls")
cat(paste(c(sampleSize,2,1),sep='\t'))
cat("\n")
cat(paste(c('#','diseased','healthy'),sep='\t'))
cat("\n")
cat(paste(cellstobeconsidered.df$group,sep='\t'))
sink()

logNormData <- seuObj@assays$SCT@data[,cellstobeconsidered]

## prep ranked Gene List
HEALTHY.Group.cell <- subset(cellstobeconsidered.df,group=='healthy')$cell
DISEASED.Group.cell <- subset(cellstobeconsidered.df,group=='diseased')$cell
HEALTHY.Group.cell.mean.expression <- apply(logNormData[,HEALTHY.Group.cell],1,mean)
DISEASED.Group.cell.mean.expression <- apply(logNormData[,DISEASED.Group.cell],1,mean)
HEALTHY.Group.cell.sd.expression <- apply(logNormData[,HEALTHY.Group.cell],1,sd)
DISEASED.Group.cell.sd.expression <- apply(logNormData[,DISEASED.Group.cell],1,sd)
signal2noise <- (DISEASED.Group.cell.mean.expression - HEALTHY.Group.cell.mean.expression)/(DISEASED.Group.cell.sd.expression + HEALTHY.Group.cell.sd.expression)
signal2noise.noNA <- signal2noise[!is.na(signal2noise)]
signal2noise.noNA <- sort(signal2noise.noNA,decreasing=T)
ranked.df <- data.frame(names(signal2noise.noNA),signal2noise.noNA)
write.table(ranked.df,file='rankedGeneList.rnk',sep='\t',row.names=F,col.names=F,quote=F)
