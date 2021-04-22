library(scater)
library(scran)


human.counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
# human.counts <- read.table("output/star/Pt0/human_genes.tsv", sep="\t", header=TRUE, row.names=1)
spikein.counts <- read.table(snakemake@input[[2]], sep="\t", header=TRUE, row.names=1)
# spikein.counts <- read.table("output/star/Pt0/spike_ins.tsv", sep="\t", header=TRUE, row.names=1)
column.names <- intersect(unique(colnames(human.counts)), unique(colnames(spikein.counts)))
# remove samples with no spike-ins

# microbe.counts <- read.table("output/Pt0/S0/D23580-cell-by-gene-count-table.tsv", sep="\t", header=TRUE, row.names=1)
microbe.counts <- read.table(snakemake@input[[3]], sep="\t", header=TRUE, row.names=1)
microbe.counts <- microbe.counts[rowSums(microbe.counts) > 0,]

column.names <- intersect(colnames(microbe.counts), column.names)
spikein.counts <- spikein.counts[, column.names]
human.counts <- human.counts[, column.names]
microbe.counts <- microbe.counts[, column.names]
# row.names=2 because 2nd row is cell
pdata <- read.table(snakemake@input[[4]], sep="\t", header = TRUE, row.names=2)
#pdata <- read.table("data/units.tsv", sep="\t", header = TRUE, row.names=2)
pdata <- pdata[column.names, ]
row.names(pdata) <- gsub("-", ".", row.names(pdata))

# create SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(counts = as.matrix(human.counts)), colData=pdata)
spike_se <- SummarizedExperiment(list(counts=spikein.counts))
altExp(sce, "spike") <- spike_se
microbe_se <- SummarizedExperiment(list(counts=microbe.counts))
altExp(sce, "microbe") <- microbe_se
sce <- computeSpikeFactors(sce, "spike")
sce <- logNormCounts(sce)
altExp(sce, "microbe") <- logNormCounts(altExp(sce, "microbe"), sizeFactors(sce))

celltype.col <- snakemake@wildcards[["celltype"]]
# celltype.col <- "status"
celltype.of.interest <- snakemake@wildcards[["celltype_of_interest"]]
# celltype.of.interest <- "D23580-infected"
celltype.comparison <- snakemake@wildcards[["celltype_comparison"]]
#celltype.comparison <- "D23580-bystander"

# now subset the cells by celltypes of interest if applicable
if (celltype.comparison != "all"){
  sce <- sce[, sce[[celltype.col]] %in% c(celltype.of.interest, celltype.comparison)]
}

min.exp <- 0
min.proportion <- .25

# get CPM using edgeR
# SingleCellExperiment::cpm(sce) <- edgeR::cpm(counts(sce))
#
# get the proportion of cells with cpm > for all possible combinations of cell group and gene
propOverXbyGroup <- apply(
  assays(altExp(sce, "microbe"))$counts,
  1,
  function(x){
  tapply(x, sce[[celltype.col]], function(x){sum(x > min.exp) / length(x)})
})

propOverXbyGroup <- reshape2::melt(
  propOverXbyGroup, varnames = c("group", "feature"), value.name = "Proportion")

# boolean for whether the proportion is > .25
propOverXbyGroup$PropOverXpct <- (propOverXbyGroup$Proportion > min.proportion)

# now we have a data.frame with four columns - group, feature, proportion and PropOver25pct

groupsOverXpct <- as.data.frame(with(
    propOverXbyGroup,
    tapply(PropOverXpct, feature, function(x){sum(x)}))
)
colnames(groupsOverXpct) <- "Groups"
groupsOverXpct$feature <- rownames(groupsOverXpct)


detectedFeatures <- names(which(tapply(
  propOverXbyGroup$Proportion, propOverXbyGroup$feature,
  function(x){sum(x > min.proportion) > 0})))

keepEndogenous <- (rownames(altExp(sce, "microbe")) %in% detectedFeatures)
altExp(sce, "microbe") <- altExp(sce, "microbe")[keepEndogenous]


gene_of_interest = snakemake@wildcards[["gene_of_interest"]]

gene.vector <- assays(sce)$logcounts[gene_of_interest,]
gene.matrix <- t(as.matrix(gene.vector))
rownames(gene.matrix) <- gene_of_interest
mat <- rbind(gene.matrix, assays(altExp(sce, "microbe"))$logcounts)

corr.res <- correlatePairs(x=mat, block=sce$plate, pairings=list(c(gene_of_interest), rownames(assays(altExp(sce, "microbe"))$logcounts)))


# out.corr <- apply(
#   assays(altExp(sce, "microbe"))$logcounts,
#   1,
#   function(x){
#   htest.obj <- cor.test(x, logcounts(sce)[gene_of_interest,], method="spearman")
#   c(htest.obj$estimate[[1]], htest.obj$p.value)
# })
#
# out.corr <- as.data.frame(t(out.corr))
#
# colnames(out.corr) <- c("rho", "pvalue")
# out.corr$FDR <- p.adjust(out.corr$pvalue)
#
# # out.corr <- cor(t(assays(altExp(sce, "microbe"))$logcounts), logcounts(sce)[gene_of_interest,], method="spearman")
# # out.corr <- out.corr[!is.na(out.corr), ]
# out.corr <- out.corr[order(out.corr$pvalue),]

write.table(corr.res, snakemake@output[[1]], sep="\t")
