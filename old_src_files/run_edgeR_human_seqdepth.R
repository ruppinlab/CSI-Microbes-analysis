library(edgeR)

counts <- read.table(snakemake@input[[1]], sep="\t", header = TRUE, row.names=1)
# counts <- read.table("output/genus_microbe_reads.tsv", sep="\t", header=TRUE, row.names=1)
pdata <- read.table(snakemake@input[[2]], sep="\t", header = TRUE, row.names=1)
pdata$plate <- as.factor(pdata$plate)
human.reads <- read.table(snakemake@input[[3]], sep="\t", header = TRUE, row.names=1)
human.reads <- human.reads[, colnames(counts)]
seq.depth <- colSums(human.reads)
stopifnot(identical(row.names(pdata), colnames(counts)))
stopifnot(identical(row.names(pdata), names(seq.depth)))
if (nlevels(pdata) > 1) {
  formula.str <- paste("~plate", snakemake@wildcards[["celltype"]], sep="+")
} else {
  formula.str <- paste0("~", snakemake@wildcards[["celltype"]])
}

formula <- as.formula(formula.str)
design <- model.matrix(formula, data=pdata)

dge <- DGEList(counts = counts, lib.size = seq.depth)
# keep all microbes with at least 10 CPM in at least .7 of the smallest group
keep <- filterByExpr(dge, design, lib.size=seq.depth)
dge <- dge[keep, , keep.lib.sizes=TRUE]
# dge <- calcNormFactors(dge)
dge <- estimateDisp(dge, design, robust=TRUE)
fit <- glmQLFit(dge, design, robust=TRUE)
glt <- glmTreat(fit, coef=ncol(design), lfc=1)
# plotMD(glt)
# abline(h=c(-1, 1), col="blue")
results <- as.data.frame(
    topTags(glt, n=Inf, adjust.method="BH", sort.by="PValue")
)
# print(results)
# print(ncol(design))
# print(glt)
write.table(results, file=snakemake@output[[1]], sep="\t", quote=FALSE)
