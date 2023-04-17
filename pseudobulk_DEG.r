celltype = "lEC"

S1.cell <- row.names(subset(meta.data.table,orig.ident=='AD1565' & cellType==celltype))
S2.cell <- row.names(subset(meta.data.table,orig.ident=='AD1961' & cellType==celltype))
S3.cell <- row.names(subset(meta.data.table,orig.ident=='AD2069' & cellType==celltype))
S4.cell <- row.names(subset(meta.data.table,orig.ident=='AD2099' & cellType==celltype))
S5.cell <- row.names(subset(meta.data.table,orig.ident=='YZ18' & cellType==celltype))
S6.cell <- row.names(subset(meta.data.table,orig.ident=='YZ19' & cellType==celltype))
S7.cell <- row.names(subset(meta.data.table,orig.ident=='YZ30' & cellType==celltype))
S8.cell <- row.names(subset(meta.data.table,orig.ident=='YZ7' & cellType==celltype))


S1.counts <- counts[,S1.cell]
S2.counts <- counts[,S2.cell]
S3.counts <- counts[,S3.cell]
S4.counts <- counts[,S4.cell]
S5.counts <- counts[,S5.cell]
S6.counts <- counts[,S6.cell]
S7.counts <- counts[,S7.cell]
S8.counts <- counts[,S8.cell]


library(Matrix)
S1.data <- rowSums(S1.counts)
S2.data <- rowSums(S2.counts)
S3.data <- rowSums(S3.counts)
S4.data <- rowSums(S4.counts)
S5.data <- rowSums(S5.counts)
S6.data <- rowSums(S6.counts)
S7.data <- rowSums(S7.counts)
S8.data <- rowSums(S8.counts)



df <- data.frame(AD1565=S1.data,AD1961=S2.data,AD2069=S3.data,AD2099=S4.data,YZ18=S5.data,
YZ19=S6.data,YZ30=S7.data,YZ7=S8.data)

coldata <- data.frame(condition = c(rep('CASE',4),rep('CTRL',4)), row.names=colnames(df))
coldata$sex <- c('female','female','male','male','female','male','female','male')
library("BiocParallel")
register(MulticoreParam(10))

#By adding variables to the design, one can control for additional variation in the counts. For example, if the condition samples are balanced across experimental batches, by including the batch factor to the design, one can increase the sensitivity for finding differences due to condition.

dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = coldata,
                              design= ~ sex + condition)

#Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# set the reference level
dds$condition <- factor(dds$condition, levels = c("CTRL","CASE"))

dds <- DESeq(dds)
res <- results(dds, name="condition_CASE_vs_CTRL", alpha=0.05)
#res <- lfcShrink(dds, coef="condition_CASE_vs_CTRL", type="ashr")
resOrdered <- res[order(res$pvalue),]

summary(res)
sum(res$padj < 0.05, na.rm=TRUE)
res.df <- data.frame(res)
res.df$gene <- row.names(res.df)
res.sig.df <- subset(res.df,padj < 0.05)

write.table(res.df,file=paste0(celltype,'.pseudobulk_results.tsv'),sep= '\t',row.names=F,col.names=T,quote=F)
write.table(res.sig.df,file=paste0(celltype,'.pseudobulk_results.padj0.05.tsv'),sep= '\t',row.names=F,col.names=T,quote=F)

geneType <- read.table(file='/work/project/xuanyu/resource/10XGenomics/refdata-gex-GRCh38-2020-A/genes.gtf.geneType.tsv',sep='\t',header=F)

protein_coding_genes <- subset(geneType, V2 == 'protein_coding')$V1

res.df.prot <- res.df[row.names(res.df) %in% protein_coding_genes,]
write.table(res.df.prot,file=paste0(celltype,'.pseudobulk_results_prot.tsv'),sep= '\t',row.names=F,col.names=T,quote=F)
res.sig.df.prot <- subset(res.df.prot,padj < 0.05)
write.table(res.sig.df.prot,file=paste0(celltype,'.pseudobulk_results_prot.padj0.05.tsv'),sep= '\t',row.names=F,col.names=T,quote=F)
save(res,file=paste0(celltype,'.res.Rdata'))
