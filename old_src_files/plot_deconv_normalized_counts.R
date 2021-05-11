library(scater)
library(scran)
library(ggplot2)

counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
# counts <- read.table("output/Pt0_genus_PathSeq_microbe_reads.tsv", sep="\t", header=TRUE, row.names=1)
pdata <- read.table(snakemake@input[[2]], sep="\t", header = TRUE, row.names=1)
row.names(pdata) <- gsub("-", ".", row.names(pdata))
# pdata <- read.table("output/Pt0_PathSeq_metadata.tsv", sep="\t", header = TRUE, row.names=1)
# filter out any celltypes with < 5 cells
celltype.col <- snakemake@wildcards[["celltype"]]
agg <- aggregate(row.names(pdata), by=list(pdata[[celltype.col]]), FUN=length)
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
#counts <- counts[rowSums(counts > min.reads) > min.cells,]
# pdata <- pdata[pdata["selection"] == "Astrocytes(HEPACAM)", ]
#pdata <- pdata[pdata["depletion_batch"] != "depleted_yes", ]
celltype <- pdata[[celltype.col]]
sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData=pdata)
sce <- computeSumFactors(sce, clusters=celltype, positive=TRUE)
sce <- logNormCounts(sce)

plotExpression(sce, x=snakemake@wildcards[["celltype"]],
    features=snakemake@wildcards[["microbe"]]) +
    theme(axis.text.x = element_text(angle = 90))

ggsave(snakemake@output[[1]])
