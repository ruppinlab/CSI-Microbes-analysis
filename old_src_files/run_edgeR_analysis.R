library(edgeR)

counts <- read.table(snakemake@input[[1]], sep="\t", header = TRUE, row.names=1)
# counts <- read.table("output/genus_microbe_reads.tsv", sep="\t", header=TRUE, row.names=1)
pdata <- read.table(snakemake@input[[2]], sep="\t", header = TRUE, row.names=1)
# filter out unknown celltypes
pdata <- pdata[pdata[snakemake@wildcards[["celltype"]]] != "unknown", ]
pdata <- droplevels(pdata)
counts <- counts[, row.names(pdata)]
print(head(pdata))
pdata$plate <- as.factor(pdata$plate)
pdata$batch <- as.factor(pdata$batch)
stopifnot(identical(row.names(pdata), colnames(counts)))
#print(snakemake@params[["formula"]])
#print(head(pdata))
formula.list <- c()
# hack to avoid "Design matrix not of full rank error"
if (nlevels(pdata$plate) >= nlevels(pdata$batch)){
  if (nlevels(pdata$plate) > 1) {
    formula.list <- append(formula.list, "plate")
  }
} else if (nlevels(pdata$batch) > 1){
  formula.list <- append(formula.list, "batch")
}

formula.list <- append(formula.list, snakemake@wildcards[["celltype"]])
formula.str <- paste0("~", paste0(formula.list, collapse="+"))


formula <- as.formula(formula.str)
design <- model.matrix(formula, data=pdata)

dge <- DGEList(counts = counts)
# keep all microbes with at least 10 CPM in at least .7 of the smallest group
keep <- filterByExpr(dge, design)
dge <- dge[keep, , keep.lib.sizes=TRUE]
dge <- calcNormFactors(dge, method=snakemake@wildcards[["norm_method"]])
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
