library(scater)
library(scran)

source(snakemake@params[["spike_functions"]])

sce <- generate_sce(microbe.file=snakemake@input[[1]], spikein.file=snakemake@input[[2]],
  pdata.file=snakemake@input[[3]], celltype.col=snakemake@wildcards[["celltype"]])

celltype.col <- snakemake@wildcards[["celltype"]]
celltype.of.interest <- snakemake@wildcards[["celltype_of_interest"]]
celltype.comparison <- snakemake@wildcards[["celltype_comparison"]]

# now subset the cells by celltypes of interest if applicable
if (celltype.comparison != "all"){
  sce <- sce[, sce[[celltype.col]] %in% c(celltype.of.interest, celltype.comparison)]
}

# filter so we are only comparing highly expressed OTUs to reduce FDR penalty
agg <- aggregate(colnames(sce), by=list(sce[[celltype.col]]), FUN=length)
min.cells <- min(agg$x)*.5
print(min.cells)
sce <- sce[rowSums(logcounts(sce) > 2) > min.cells, ]
# celltype.of.interest <- snakemake@wildcards[["celltype_of_interest"]]
groups <- sce[[celltype.col]]
lfc <- as.numeric(snakemake@wildcards[["lfc"]])
pval.type <- snakemake@wildcards[["pval.type"]]
block <- sce[[snakemake@wildcards[["block"]]]]
# calculate markers using t-test
t.test.markers <- findMarkers(sce, groups=groups, lfc=lfc, pval.type=pval.type, block=block)
t.test.df <- t.test.markers[[celltype.of.interest]]
write.table(t.test.df, file=snakemake@output[[1]], sep="\t")
# calculate markers using wilcoxon rank sum test
wilcox.markers <- findMarkers(sce, groups=groups, test="wilcox", lfc=lfc, pval.type=pval.type, block=block)
wilcox.df <- wilcox.markers[[celltype.of.interest]]
write.table(wilcox.df, file=snakemake@output[[2]], sep="\t")
