library(scran)
library(scater)

generate_sce <- function(gene.file, spike.file, pdata.file="data/units.tsv") {
  human.counts <- read.table(gene.file, sep="\t", header=TRUE, row.names=1)
  spikein.counts <- read.table(spike.file, sep="\t", header=TRUE, row.names=1)
  pdata <- read.table(pdata.file, sep="\t", header = TRUE, row.names=3)
  row.names(pdata) <- gsub("-", ".", row.names(pdata))
  column.names <- intersect(unique(colnames(human.counts)), unique(rownames(pdata)))
  column.names <- column.names[!(column.names %in% c("B1_B003648", "P5_B000420"))]
  spikein.counts <- spikein.counts[, column.names]
  human.counts <- human.counts[, column.names]
  pdata <- pdata[column.names, ]
  sce <- SingleCellExperiment(assays = list(counts = as.matrix(human.counts)), colData=pdata)
  spike_se <- SummarizedExperiment(list(counts=spikein.counts))
  altExp(sce, "spike") <- spike_se
  return(sce)
}


sce.list <- mapply(generate_sce, snakemake@input[["genes"]], snakemake@input[["spikes"]])
sce <- do.call(cbind, sce.list)
sce <- computeSpikeFactors(sce, "spike")
sce <- logNormCounts(sce)

# focus only on tumor cells
sce <- sce[,sce$Tumor == "Tumor"]

saveRDS(sce, snakemake@output[[1]])
