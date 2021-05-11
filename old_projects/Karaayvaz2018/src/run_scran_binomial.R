library(dplyr)
library(scater)
library(scran)
library(ggplot2)


counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
# counts <- read.table("output/Pt0_genus_PathSeq_microbe_reads.tsv", sep="\t", header=TRUE, row.names=1)
pdata <- read.table(snakemake@input[[2]], sep="\t", header = TRUE, row.names=1)
row.names(pdata) <- gsub("-", ".", row.names(pdata))
# pdata <- read.table("output/Pt0_PathSeq_metadata.tsv", sep="\t", header = TRUE, row.names=1)
# filter out any celltypes with < 5 cells
celltype.col <- "QC_status"
celltype.of.interest <- "Passed"
agg <- aggregate(row.names(pdata), by=list(pdata[[celltype.col]]), FUN=length)
min.cells <- min(agg$x)*.7
sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData=pdata)
sce <- sce[, colSums(counts(sce)) > 0]
#sce <- computeSumFactors(sce, cluster=celltype)
sce <- logNormCounts(sce, size_factors=librarySizeFactors(sce))
# filter so we are only comparing highly expressed OTUs to reduce FDR penalty
sce <- sce[rowSums(logcounts(sce) > 2) > min.cells, ]
# cluster by celltype
# if (nlevels(pdata$plate) >= nlevels(pdata$batch)){
#   block <- pdata$plate
# } else {
#   block <- pdata$batch
# }
groups <- colData(sce)[[celltype.col]]
lfc <- 1
pval.type <- "all"

# t <- findMarkers(sce, groups=groups, lfc=lfc, pval.type=pval.type, block=block)
# df <- t[[celltype.of.interest]]
# write.table(df, file=snakemake@output[[1]], sep="\t")
# wilcox <- findMarkers(sce, test="wilcox", groups=groups, lfc=lfc, pval.type=pval.type, block=block)
# df <- wilcox[[celltype.of.interest]]
# write.table(df, file=snakemake@output[[2]], sep="\t")
binom <- findMarkers(sce, test="binom", groups=groups, lfc=lfc, pval.type=pval.type)
df <- binom[[celltype.of.interest]]
write.table(df, file=snakemake@output[[1]], sep="\t")
# combined <- multiMarkerStats(
#     t=t,
#     wilcox=wilcox,
#     binom=binom
# )
#
# df <- combined[[celltype.of.interest]]
# write.table(df, file=snakemake@output[[1]], sep="\t")
