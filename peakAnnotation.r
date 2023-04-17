gr <- rtracklayer::import('/work/project/xuanyu/resource/10XGenomics/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz')
txdb <- makeTxDbFromGRanges(gr)
annoData = genes(txdb)
geneid_genename <- as.data.frame(gr)[,c('gene_id','gene_name')]
geneid_genename <- unique(geneid_genename)
library(plyr)
genename <- mapvalues(annoData$gene_id,from=geneid_genename$gene_id,to=geneid_genename$gene_name)
annoData$gene_name <- genename
seqlevelsStyle(annoData) <- 'UCSC'

library(Signac)
library(Seurat)
load('seuObj.clean.Rdata')
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
peaks <- seuObj@assays$peaks@ranges

# Peak distribution over different genomic features.
     The distribution will be calculated by geneLevel, ExonIntron, and
     Exons The geneLevel will be categorized as promoter region, gene
     body, gene downstream and distal intergenic region. The ExonIntron
     will be categorized as exon, intron and intergenic. The Exons will
     be categorized as 5' UTR, 3'UTR and CDS. 
     The precedence will follow the order of labels defination. For example, for
     ExonIntron, if a peak overlap with both exon and intron, and exon
     is specified before intron, then only exon will be incremented for
     the same example.


out <- genomicElementDistribution(peaks, 
                           TxDb = txdb,
                           promoterRegion=c(upstream=2000, downstream=500),
                           geneDownstream=c(upstream=0, downstream=2000),
                           promoterLevel=list(
                         # from 5' -> 3', fixed precedence 3' -> 5'
                             breaks = c(-2000, -1000, -500, 0, 500),
                             labels = c("upstream 1-2Kb", "upstream 0.5-1Kb", 
                                        "upstream <500b", "TSS - 500b"),
                             colors = c("#FFE5CC", "#FFCA99", 
                                        "#FFAD65", "#FF8E32")))
#To obtain peaks with nearest bi-directional promoters
              within 5kb upstream and 3kb downstream of TSS, set output
              = "nearestBiDirectionalPromoters" and bindingRegion =
              c(-5000, 3000)
            • To obtain peaks within 5kb upstream and up to 3kb
              downstream of TSS within the gene body, set
              output="overlapping", FeatureLocForDistance="TSS" and
              bindingRegion = c(-5000, 3000)

            • To obtain peaks up to 5kb upstream within the gene body
              and 3kb downstream of gene/Exon End, set
              output="overlapping", FeatureLocForDistance="geneEnd" and
              bindingRegion = c(-5000, 3000)

            • To obtain peaks from 5kb upstream to 3kb downstream of
              genes/Exons, set output="overlapping", bindingType =
              "fullRange" and bindingRegion = c(-5000, 3000)
                                       
overlaps.anno.promoter <- annoPeaks(out$peaks, 
                                     annoData=annoData, 
                                     bindingType="startSite",
                                     bindingRegion=c(-2000, 500),
                                     ignore.peak.strand = TRUE,
                                     select = "bestOne"
                                     )
overlaps.anno.promoter <- as.data.frame(subset(overlaps.anno.promoter,geneLevel=='promoter'))
overlaps.anno.downstream <- annoPeaks(out$peaks, 
                                     annoData=annoData, 
                                     bindingType="endSite",
                                     bindingRegion=c(0, 2000),
                                     ignore.peak.strand = TRUE,
                                     select = "bestOne"
                                     )
overlaps.anno.downstream <-  as.data.frame(subset(overlaps.anno.downstream,geneLevel=='geneDownstream'))

overlaps.anno.genebody <- annoPeaks(out$peaks, 
                                     annoData=annoData, 
                                     bindingType="fullRange",
                                     bindingRegion=c(0, 1),
                                     ignore.peak.strand = TRUE,
                                     select = "bestOne"
                                     )
overlaps.anno.genebody <-  as.data.frame(subset(overlaps.anno.genebody,geneLevel=='geneBody'))

# we added upstream annotation
overlaps.anno.upstream <- annoPeaks(out$peaks, 
                                     annoData=annoData, 
                                    bindingType="startSite",
                                     bindingRegion=c(-5000, 1),
                                     ignore.peak.strand = TRUE,
                                     select = "bestOne"
                                     )
overlaps.anno.upstream <-  as.data.frame(subset(overlaps.anno.upstream,geneLevel=='distalIntergenic'))
overlaps.anno.upstream$geneLevel <- gsub('distalIntergenic','upstream',overlaps.anno.upstream$geneLevel)
overlaps.anno.upstream.id <- paste0(paste0(overlaps.anno.upstream$seqnames,':',overlaps.anno.upstream$start),'-',overlaps.anno.upstream$end)


overlaps.anno.distalIntergenic <- as.data.frame(subset(out$peaks,geneLevel=='distalIntergenic'))
overlaps.anno.distalIntergenic.id <- paste0(paste0(overlaps.anno.distalIntergenic$seqnames,':',overlaps.anno.distalIntergenic$start),'-',overlaps.anno.distalIntergenic$end)

overlaps.anno.distalIntergenic <- overlaps.anno.distalIntergenic[!overlaps.anno.distalIntergenic.id %in% overlaps.anno.upstream.id,]

addcol <- c('peak','feature','feature.ranges.start','feature.ranges.end','feature.ranges.width','feature.strand','distance','insideFeature','distanceToSite','gene_id','gene_name')
newMatrix <- matrix(rep('-',nrow(overlaps.anno.distalIntergenic)*length(addcol)),nrow=nrow(overlaps.anno.distalIntergenic), ncol=length(addcol))
newDataf <- as.data.frame(newMatrix)
colnames(newDataf) <- addcol
overlaps.anno.distalIntergenic <- cbind(overlaps.anno.distalIntergenic,newDataf)

anno.combined <- rbind(overlaps.anno.promoter,overlaps.anno.downstream,overlaps.anno.genebody,overlaps.anno.upstream,overlaps.anno.distalIntergenic)
anno.combined <- anno.combined[with(anno.combined, order(seqnames, start, end)),]
anno.combined.id <- paste0(paste0(anno.combined$seqnames,':',anno.combined$start),'-',anno.combined$end)
anno.combined$peak_id <- anno.combined.id
peaks_df <- as.data.frame(out$peaks)
peaks_df.id <- paste0(paste0(peaks_df$seqnames,':',peaks_df$start),'-',peaks_df$end)
peaks_df$peak_id <- peaks_df.id
unAnnotPeaks <-  peaks_df[!peaks_df$peak_id  %in% anno.combined$peak_id,]
print(paste0("The number of not annotated peaks: ",nrow(unAnnotPeaks)))
newMatrix2 <- matrix(rep('-',nrow(unAnnotPeaks)*length(addcol)),nrow=nrow(unAnnotPeaks), ncol=length(addcol))
newDataf2 <- as.data.frame(newMatrix2)
colnames(newDataf2) <- addcol
unAnnotPeaks <- cbind(unAnnotPeaks[,1:9],newDataf2)
unAnnotPeaks.id <- paste0(paste0('unannot_',unAnnotPeaks$seqnames,':',unAnnotPeaks$start),'-',unAnnotPeaks$end)
unAnnotPeaks$peak_id <- unAnnotPeaks.id

anno.combined.all <- rbind(anno.combined,unAnnotPeaks)
anno.combined.all <- anno.combined.all[with(anno.combined.all, order(seqnames, start, end)),]
write.table(anno.combined.all,col.names=T,row.names=F,sep='\t',file='all_annotated_peaks.tsv',quote=F)