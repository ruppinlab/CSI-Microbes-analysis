library(dplyr)
library(scater)
library(scran)
library(ggplot2)


counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
# counts <- read.table("output/Pt0_genus_PathSeq_microbe_reads.tsv", sep="\t", header=TRUE, row.names=1)
pdata <- read.table(snakemake@input[[2]], sep="\t", header = TRUE, row.names=1)
human.reads <- read.table(snakemake@input[[3]], sep="\t", header = TRUE, row.names=1)
# human.reads <- read.table("output/star-readcounts.tsv", sep="\t", header = TRUE, row.names=1)


# we want to subset the data so that we get a subset of the samples
infected.plate <- snakemake@wildcards[["infected_plate"]]
control.plate <- snakemake@wildcards[["control_plate"]]
pdata <- pdata[((pdata$infected == "infected") & (pdata$plate == infected.plate)) | ((pdata$infected == "uninfected") & (pdata$plate == control.plate)), ]

# filter out any celltypes with < 5 cells
celltype.col <- snakemake@wildcards[["celltype"]]
celltype.of.interest <- snakemake@wildcards[["celltype_of_interest"]]
agg <- aggregate(row.names(pdata), by=list(pdata[[celltype.col]]), FUN=length)
celltypes.to.keep <- agg[agg$x >= 5,]$Group.1
pdata <- pdata[pdata[[celltype.col]] %in% celltypes.to.keep, ]
# filter out unknown celltypes
pdata <- pdata[pdata[[celltype.col]] != "unknown", ]
pdata <- droplevels(pdata)
counts <- counts[, row.names(pdata)]
# filter out any OTUs without the minimum number of reads in .7 of the smallest group
# before normalization -> reduce sparsity to manageable amount
min.reads <- 2
min.cells <- min(agg$x)*.7
counts <- counts[rowSums(counts > min.reads) > min.cells,]

human.reads <- human.reads[, colnames(counts)]
spike.prefix <- snakemake@params[["spike"]]
spike.reads <- human.reads[grep(spike.prefix, row.names(human.reads)),]

sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData=pdata)
spike_se <- SummarizedExperiment(list(counts=spike.reads))
altExp(sce, "spike") <- spike_se
sce <- computeSpikeFactors(sce, "spike")
# pdata <- pdata[pdata["selection"] == "Astrocytes(HEPACAM)", ]
#pdata <- pdata[pdata["depletion_batch"] != "depleted_yes", ]
celltype <- pdata[[celltype.col]]
sce <- logNormCounts(sce)
# filter so we are only comparing highly expressed OTUs to reduce FDR penalty
sce <- sce[rowSums(logcounts(sce) > 2) > min.cells, ]
# cluster by celltype
# if (nlevels(pdata$plate) >= nlevels(pdata$batch)){
#   block <- pdata$plate
# } else {
#   block <- pdata$batch
# }
groups <- celltype
lfc <- 1
pval.type <- "some"
t <- findMarkers(sce, groups=groups, lfc=lfc, pval.type=pval.type)
df <- t[[celltype.of.interest]]
write.table(df, file=snakemake@output[[1]], sep="\t")
wilcox <- findMarkers(sce, test="wilcox", groups=groups, lfc=lfc, pval.type=pval.type)
df <- wilcox[[celltype.of.interest]]
write.table(df, file=snakemake@output[[2]], sep="\t")
# combined <- multiMarkerStats(
#     t=t,
#     wilcox=wilcox,
#     binom=binom
# )
#
# df <- combined[[celltype.of.interest]]
# write.table(df, file=snakemake@output[[1]], sep="\t")
