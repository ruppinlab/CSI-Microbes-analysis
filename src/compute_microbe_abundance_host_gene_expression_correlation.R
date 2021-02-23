library(scater)
library(scran)


human.counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
# human.counts <- read.table("output/star/Pt0/human_genes.tsv", sep="\t", header=TRUE, row.names=1)
spikein.counts <- read.table(snakemake@input[[2]], sep="\t", header=TRUE, row.names=1)
# spikein.counts <- read.table("output/star/Pt0/spike_ins.tsv", sep="\t", header=TRUE, row.names=1)
microbe.counts <- read.table(snakemake@input[[3]], sep="\t", header=TRUE, row.names=1)
# microbe.counts <- read.table("output/Pt0/genus_PathSeq_Bacteria_reads.tsv", sep="\t", header=TRUE, row.names=1)
pdata <- read.table(snakemake@input[[4]], sep="\t", header = TRUE, row.names=1)
# pdata <- read.table("data/units.tsv", sep="\t", header = TRUE, row.names=2)
row.names(pdata) <- gsub("-", ".", row.names(pdata))
# sometimes there are cells in the human.counts and spikein.counts that are
# not in microbe.counts and pdata because they were filtered by the original authors
# also remove any cells with 0 spike-in reads - c("B1_B003648", "P5_B000420")
# TODO automatically remove cells with 0 spike-in reads
column.names <- intersect(unique(colnames(human.counts)), unique(colnames(microbe.counts)))
column.names <- column.names[!(column.names %in% c("B1_B003648", "P5_B000420"))]
#print(dim(column.names))
spikein.counts <- spikein.counts[, column.names]
human.counts <- human.counts[, column.names]
microbe.counts <- microbe.counts[, column.names]
pdata <- pdata[column.names, ]

# get table to convert from tax_id to microbe name
# tax.map <- read.table("output/Pt0/tax_id_map_Bacteria_PathSeq.tsv", sep="\t", header=TRUE)
tax.map <- read.table(snakemake@input[[5]], sep="\t", header=TRUE)
# create SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(counts = as.matrix(human.counts)), colData=pdata)
spike_se <- SummarizedExperiment(list(counts=spikein.counts))
altExp(sce, "spike") <- spike_se
microbe_se <- SummarizedExperiment(list(counts=microbe.counts))
altExp(sce, "microbe") <- microbe_se
sce <- computeSpikeFactors(sce, "spike")
sce <- logNormCounts(sce)
altExp(sce, "microbe") <- logNormCounts(altExp(sce, "microbe"), sizeFactors(sce))

rownames(altExp(sce, "microbe")) <- lapply(rownames(altExp(sce, "microbe")), function(x) tax.map[tax.map$tax_id == x, "name"])

celltype.col <- snakemake@wildcards[["celltype"]]
#celltype.col <- "infected"
celltype.of.interest <- snakemake@wildcards[["celltype_of_interest"]]
#celltype.of.interest <- "infected"
celltype.comparison <- snakemake@wildcards[["celltype_comparison"]]
#celltype.comparison <- "bystander"

# now subset the cells by celltypes of interest if applicable
if (celltype.comparison != "all"){
  sce <- sce[, sce[[celltype.col]] %in% c(celltype.of.interest, celltype.comparison)]
}

#sce$group <- sce$status

# get CPM using edgeR
SingleCellExperiment::cpm(sce) <- edgeR::cpm(counts(sce))

# if there is only one unique value in sce[[celltype.col]], this function returns a list
# if there is > 1 unique value in sce[[celltype.col]], this function returns a matrix
# get the proportion of cells with cpm > for all possible combinations of cell group and gene
propOver10byGroup <- apply(
  cpm(sce),
  1,
  function(x){
  tapply(x, sce[[celltype.col]], function(x){sum(x > 10) / length(x)})
})

if (length(unique(sce[[celltype.col]])) == 1){
  propOver10byGroup <- t(as.matrix(propOver10byGroup))
  rownames(propOver10byGroup) <- unique(sce[[celltype.col]])
}
#print(propOver10byGroup)
propOver10byGroup <- reshape2::melt(
  propOver10byGroup, varnames = c("group", "feature"), value.name = "Proportion")
#print(propOver10byGroup)
# boolean for whether the proportion is > .25
propOver10byGroup$PropOver25pct <- (propOver10byGroup$Proportion > 0.25)

# now we have a data.frame with four columns - group, feature, proportion and PropOver25pct

groupsOver25pct <- as.data.frame(with(
    propOver10byGroup,
    tapply(PropOver25pct, feature, function(x){sum(x)}))
)
colnames(groupsOver25pct) <- "Groups"
groupsOver25pct$feature <- rownames(groupsOver25pct)

detectedFeatures <- names(which(tapply(
  propOver10byGroup$Proportion, propOver10byGroup$feature,
  function(x){sum(x > 0.25) > 0})))

keepEndogenous <- (rownames(sce) %in% detectedFeatures)
sce <- sce[keepEndogenous]

microbe <- snakemake@wildcards[["microbe"]]

microbe.vector <- assays(altExp(sce, "microbe"))$logcounts[microbe,]
microbe.matrix <- t(as.matrix(microbe.vector))
rownames(microbe.matrix) <- microbe
mat <- rbind(microbe.matrix, assays(sce)$logcounts)

corr.res <- correlatePairs(x=mat, block=sce$plate, pairings=list(c(microbe), rownames(assays(sce)$logcounts)))

# for each human gene X, compute the correlation between X and the microbe_of_interest for each plate
# out.corr <- apply(
#   assays(sce)$logcounts,
#   1,
#   function(x){
#     out <- rbind(x, assays(altExp(sce, "microbe"))$logcounts[microbe,])
#     #rownames(out) <- c("x", "microbe")
#     out <- vector("list", length(unique(sce$plate)))
#     pval <- vector("list", length(unique(sce$plate)))
#     i <- 1
#     for (p in unique(sce$plate)){
#       htest.obj <- cor.test(x[sce$plate == p], assays(altExp(sce, "microbe"))$logcounts[microbe, sce$plate == p], method="spearman")
#       #print(c(htest.obj$estimate[[1]], htest.obj$p.value))
#       out[[i]] <- htest.obj$estimate[[1]]
#       pval[[i]] <- htest.obj$p.value
#       i <- i+1
#     }
#     return(c(median(unlist(out)), mean(unlist(pval))))
#   }
# )

# htest.obj <- cor.test(x, assays(altExp(sce, "microbe"))$logcounts[microbe,], method="spearman")
# c(htest.obj$estimate[[1]], htest.obj$p.value)
# # out.corr <- apply(
# #   assays(altExp(sce, "microbe"))$logcounts,
# #   1,
# #   function(x){
# #   htest.obj <- cor.test(x, logcounts(sce)[gene_of_interest,], method="spearman")
# #   c(htest.obj$estimate[[1]], htest.obj$p.value)
# # })
#
# # htest.obj <- cor.test(, logcounts(sce)[gene_of_interest,], method="spearman")
# # c(htest.obj$estimate[[1]], htest.obj$p.value)
# # htest.obj <- cor.test(x, logcounts(sce)[gene_of_interest,], method="spearman")
# # c(htest.obj$estimate[[1]], htest.obj$p.value)
#
#
# out.corr <- cor(t(logcounts(sce)), assays(altExp(sce, "microbe"))$logcounts[microbe,], method="spearman")
# out.corr <- t(out.corr)
# colnames(out.corr) <- c("rho", "pval")
# out.corr <- out.corr[order(out.corr[,"rho"]),]
write.table(corr.res, snakemake@output[[1]])
