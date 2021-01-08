
# generate a SingleCellExperiment object
generate_sce <- function(microbe.file, spikein.file, pdata.file, celltype.col) {
  microbe.counts <- read.table(microbe.file, sep="\t", header=TRUE, row.names=1)
  spikein.counts <- read.table(spikein.file, sep="\t", header=TRUE, row.names=1)
  column.names <- intersect(unique(colnames(microbe.counts)), unique(colnames(spikein.counts)))
  # remove samples with no spike-ins
  column.names <- column.names[!(column.names %in% c("B1_B003648", "P5_B000420"))]
  spikein.counts <- spikein.counts[, column.names]
  microbe.counts <- microbe.counts[, column.names]
  counts <- microbe.counts
  pdata <- read.table(pdata.file, sep="\t", header = TRUE, row.names=1)
  pdata <- pdata[column.names, ]
  row.names(pdata) <- gsub("-", ".", row.names(pdata))
  #print(pdata)
  #print(celltype.col)
  #print(pdata[[celltype.col]])
  agg <- aggregate(row.names(pdata), by=list(pdata[[celltype.col]]), FUN=length)
  # filter out any celltypes with < 5 cells
  celltypes.to.keep <- agg[agg$x >= 5,]$Group.1
  pdata <- pdata[pdata[[celltype.col]] %in% celltypes.to.keep, ]
  # filter out unknown celltypes
  pdata <- pdata[pdata[[celltype.col]] != "Unknown", ]
  pdata <- droplevels(pdata)
  counts <- counts[, row.names(pdata)]
  spikein.counts <- spikein.counts[, row.names(pdata)]

  #rownames(counts) <- lapply(rownames(counts), function(x) tax.map[tax.map$tax_id == x, "name"])
  sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData=pdata)
  #sce <- sce[, colSums(counts(sce)) > 0]
  celltype <- sce[[celltype.col]]
  spike_se <- SummarizedExperiment(list(counts=spikein.counts))
  altExp(sce, "spike") <- spike_se
  sce <- computeSpikeFactors(sce, "spike")
  sce <- logNormCounts(sce)
  return(sce)
}
