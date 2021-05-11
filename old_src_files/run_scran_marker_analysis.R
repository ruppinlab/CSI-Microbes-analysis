library(scater)
library(scran)
library(ggplot2)


counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
# counts <- read.table("output/Pt0_species_PathSeq_Bacteria_reads.tsv", sep="\t", header=TRUE, row.names=1)
pdata <- read.table(snakemake@input[[2]], sep="\t", header = TRUE, row.names=1)
# pdata <- read.table("output/Pt0_species_PathSeq_Bacteria_metadata.tsv", sep="\t", header = TRUE, row.names=1)
row.names(pdata) <- gsub("-", ".", row.names(pdata))

celltype.col <- snakemake@wildcards[["celltype"]]
celltype.of.interest <- snakemake@wildcards[["celltype_of_interest"]]
celltype.comparison <- snakemake@wildcards[["celltype_comparison"]]
if (celltype.comparison != "all"){
  pdata <- pdata[pdata[[celltype.col]] %in% c(celltype.of.interest, celltype.comparison), ]
  counts <- counts[, row.names(pdata)]
}

agg <- aggregate(row.names(pdata), by=list(pdata[[celltype.col]]), FUN=length)
# filter out any celltypes with < 5 cells
celltypes.to.keep <- agg[agg$x >= 5,]$Group.1
pdata <- pdata[pdata[[celltype.col]] %in% celltypes.to.keep, ]
# filter out unknown celltypes
pdata <- pdata[pdata[[celltype.col]] != "Unknown", ]
pdata <- droplevels(pdata)
counts <- counts[, row.names(pdata)]
# filter out any OTUs without the minimum number of reads in .7 of the smallest group
# before normalization -> reduce sparsity to manageable amount
min.reads <- 2
min.cells <- min(agg$x)*.5
counts <- counts[rowSums(counts > min.reads) > min.cells,]
sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData=pdata)
# get the cell-type for each cell
celltype <- sce[[celltype.col]]
# compute the normalization factor by clustering cells from the same cell-type together
sce <- computeSumFactors(sce, min.mean=1)
sce <- logNormCounts(sce)
# Further filter OTUs so we are only comparing highly expressed OTUs to reduce FDR penalty
sce <- sce[rowSums(logcounts(sce) > 2) > min.cells, ]
groups <- celltype
lfc <- as.numeric(snakemake@wildcards[["lfc"]])
pval.type <- snakemake@wildcards[["pvaltype"]]
block <- pdata[[snakemake@wildcards[["block"]]]]
# calculate markers using t-test
t <- findMarkers(sce, groups=groups, lfc=lfc, pval.type=pval.type, block=block)
df <- t[[celltype.of.interest]]
write.table(df, file=snakemake@output[[1]], sep="\t")
# calculate markers using Wilcoxon-rank sum test
wilcox <- findMarkers(sce, test="wilcox", groups=groups, lfc=lfc, pval.type=pval.type, block=block)
df <- wilcox[[celltype.of.interest]]
write.table(df, file=snakemake@output[[2]], sep="\t")
