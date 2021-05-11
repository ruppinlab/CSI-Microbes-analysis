library(scater)
library(scran)
library(ggplot2)


#counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
counts <- read.table("output/Pt0_species_PathSeq_Bacteria_reads.tsv", sep="\t", header=TRUE, row.names=1)
#pdata <- read.table(snakemake@input[[2]], sep="\t", header = TRUE, row.names=1)
pdata <- read.table("output/Pt0_species_PathSeq_Bacteria_metadata.tsv", sep="\t", header = TRUE, row.names=1)
row.names(pdata) <- gsub("-", ".", row.names(pdata))

#celltype.col <- snakemake@wildcards[["celltype"]]
celltype.col <- "infected"
# celltype.of.interest <- snakemake@wildcards[["celltype_of_interest"]]
celltype.of.interest <- "infected"
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
min.cells <- min(agg$x)*.7
counts <- counts[rowSums(counts > min.reads) > min.cells,]
sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData=pdata)
# get the cell-type for each cell
celltype <- sce[[celltype.col]]
# aggregateAcrossCells or sumCountsAcrossCells
# summed <- aggregateAcrossCells(sce, ids=colData(sce)[,c(celltype.col, snakemake@wildcards[["block"]])])
summed <- aggregateAcrossCells(sce, ids=colData(sce)[,c(celltype.col, "plate")])
summed <- summed[, summed$infected %in% c("bystander", "infected")]
y <- DGEList(counts(summed), samples=colData(summed))
y <- calcNormFactors(y)
design <- model.matrix(~factor(plate) + factor(infected), y$samples)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
res <- glmQLFTest(fit, coef=ncol(design))
topTags(res)
# label <- "infected"
# contains pseudo bulk
current <- summed[, label==summed[[celltype.col]]]
design <-  model.matrix(~factor(plate) + factor(infected), data=colData(summed))

# sample is like the patient - controlling for this
# we want the
de.results <- pseudoBulkDGE(summed,
    sample=summed$plate,
    label=summed[[celltype.col]],
    design=design

    # 'condition' sets the group size for filterByExpr(),
    # to perfectly mimic our previous manual analysis.
    condition=targets$tomato
    #     sample=summed.filt$sample,
)
