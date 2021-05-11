library(scater)
library(scran)
library(ggplot2)


counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
# counts <- read.table("output/Pt0_genus_PathSeq_microbe_reads.tsv", sep="\t", header=TRUE, row.names=1)
pdata <- read.table(snakemake@input[[2]], sep="\t", header = TRUE, row.names=1)
# pdata <- read.table("output/Pt0_PathSeq_metadata.tsv", sep="\t", header = TRUE, row.names=1)
row.names(pdata) <- gsub("-", ".", row.names(pdata))
# spike-in reads
spike.prefix <- snakemake@params[["spike"]]
spike.reads <- counts[grep(spike.prefix, row.names(counts)),]
# rename some wildcards
celltype.col <- snakemake@wildcards[["celltype"]]
celltype.of.interest <- snakemake@wildcards[["celltype_of_interest"]]
agg <- aggregate(row.names(pdata), by=list(pdata[[celltype.col]]), FUN=length)
# filter out any celltypes with < 5 cells
#celltypes.to.keep <- agg[agg$x >= 5,]$Group.1
#pdata <- pdata[pdata[[celltype.col]] %in% celltypes.to.keep, ]
# filter out unknown celltypes
pdata <- pdata[pdata[[celltype.col]] != "Unknown", ]
pdata <- droplevels(pdata)
counts <- counts[, row.names(pdata)]
spike.reads <- spike.reads[, row.names(pdata)]
# filter out any OTUs without the minimum number of reads in .7 of the smallest group
# before normalization -> reduce sparsity to manageable amount
min.reads <- 1
min.cells <- min(agg$x)*.1
counts <- counts[rowSums(counts > min.reads) > min.cells,]
sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData=pdata)
#sce <- sce[, colSums(counts(sce)) > 0]
celltype <- sce[[celltype.col]]
spike_se <- SummarizedExperiment(list(counts=spike.reads))
altExp(sce, "spike") <- spike_se
sce <- computeSpikeFactors(sce, "spike")
sce <- logNormCounts(sce)
print(sce)
# filter so we are only comparing highly expressed OTUs to reduce FDR penalty
#sce <- sce[rowSums(logcounts(sce) > 2) > min.cells, ]
print(sce)
groups <- celltype
lfc <- as.numeric(snakemake@wildcards[["lfc"]])
pval.type <- snakemake@wildcards[["pval.type"]]
print(snakemake@wildcards[["block"]])
block <- pdata[[snakemake@wildcards[["block"]]]]
print(block)
# calculate markers using t-test
t.test.markers <- findMarkers(sce, groups=groups, lfc=lfc, pval.type=pval.type, block=block)
print(t.test.markers)
t.test.df <- t.test.markers[[celltype.of.interest]]
print(t.test.df)
t.test.df <- t.test.df[grep(snakemake@wildcards[["kingdom"]], row.names(t.test.df)),]
t.test.df[["FDR"]] <- p.adjust(t.test.df[["p.value"]], method="fdr")
row.names(t.test.df) <- gsub(paste0(snakemake@wildcards[["kingdom"]], "-"), "", row.names(t.test.df))
write.table(t.test.df, file=snakemake@output[[1]], sep="\t")
#wilcox.markers <- findMarkers(sce, groups=celltype, test="wilcox", lfc=1, block=block)
#write.table(wilcox.markers[[celltype.of.interest]], file=snakemake@output[[2]])
wilcox.markers <- findMarkers(sce, groups=groups, test="wilcox", lfc=lfc, pval.type=pval.type, block=block)
wilcox.df <- wilcox.markers[[celltype.of.interest]]
wilcox.df <- wilcox.df[grep(snakemake@wildcards[["kingdom"]], row.names(wilcox.df)),]
wilcox.df[["FDR"]] <- p.adjust(wilcox.df[["p.value"]], method="fdr")
row.names(wilcox.df) <- gsub(paste0(snakemake@wildcards[["kingdom"]], "-"), "", row.names(wilcox.df))
write.table(wilcox.df, file=snakemake@output[[2]], sep="\t")
