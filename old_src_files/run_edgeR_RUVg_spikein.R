library(edgeR)
library(RUVSeq)


counts <- read.table(snakemake@input[[1]], sep="\t", header = TRUE, row.names=1)
# counts <- read.table("output/genus_microbe_reads.tsv", sep="\t", header=TRUE, row.names=1)
pdata <- read.table(snakemake@input[[2]], sep="\t", header = TRUE, row.names=1)
pdata$batch <- as.factor(pdata$plate)
human.reads <- read.table(snakemake@input[[3]], sep="\t", header = TRUE, row.names=1)
human.reads <- human.reads[, colnames(counts)]
spike.prefix <- snakemake@params[["spike"]]
spike.reads <- human.reads[grep(spike.prefix, row.names(human.reads)),]
spike.genes <- rownames(human.reads)[grep(spike.prefix, row.names(human.reads))]
combined.counts <- rbind(counts, spike.reads)
stopifnot(identical(row.names(pdata), colnames(counts)))
if (nlevels(pdata$batch) > 1) {
  formula.str <- paste("~batch", snakemake@wildcards[["celltype"]], sep="+")
} else {
  formula.str <- paste0("~", snakemake@wildcards[["celltype"]])
}

formula <- as.formula(formula.str)
design <- model.matrix(formula, data=pdata)

dge <- DGEList(counts = counts)
# keep all microbes with at least 10 CPM in at least .7 of the smallest group
keep <- filterByExpr(dge, design)
dge <- dge[keep, , keep.lib.sizes=TRUE]
# calculate norm factors using the spike-in sequences
#N <- colSums(counts)
#print(N)
dge <- calcNormFactors(dge, method="upperquartile")
#print(nf)
#dge$samples$norm.factors <- calcNormFactors(spike.reads, lib.size=N)
#dge$samples$norm.factors <- calcNormFactors(dge[grep(spike.prefix, row.names(dge$counts)),])$samples$norm.factors
dge <- estimateDisp(dge, design, robust=TRUE)
# given raw counts, fits negative
fit <- glmQLFit(dge, design, robust=TRUE)
glt <- glmTreat(fit, coef=ncol(design), lfc=1)
# plotMD(glt)
# abline(h=c(-1, 1), col="blue")
results <- as.data.frame(
    topTags(glt, n=Inf, adjust.method="BH", sort.by="PValue")
)

write.table(results, file=snakemake@output[[1]], sep="\t", quote=FALSE)
