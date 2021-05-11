library(metagenomeSeq)

counts <- read.table(snakemake@input[[1]], sep="\t", header = TRUE, row.names=1)
# counts <- read.table("output/genus_microbe_Pt0_reads.tsv", sep="\t", header = TRUE, row.names=1)
pdata <- read.table(snakemake@input[[2]], sep="\t", header = TRUE, row.names=1)
# pdata <- read.table("output/genus_microbe_Pt0_metadata.tsv", sep="\t", header = TRUE, row.names=1)
pdata$plate <- as.factor(pdata$plate)
stopifnot(identical(row.names(pdata), colnames(counts)))

# if (nlevels(pdata) > 1) {
#   formula.str <- paste("~plate", snakemake@params["formula"], snakemake@wildcards[["celltype"]], sep="+")
# } else {
#   formula.str <- paste0("~", snakemake@params["formula"], "+", snakemake@wildcards[["celltype"]])
# }
formula <- as.formula(formula.str)
design <- model.matrix(formula, data=pdata)

exp <- newMRexperiment(counts, phenoData = AnnotatedDataFrame(pdata))
exp <- filterData(exp, present=30, depth=1)
# getting two warnings with wrenchNorm
# 1: glm.fit: algorithm did not converge
# 2: glm.fit: fitted probabilities numerically 0 or 1 occurred
# exp <- wrenchNorm(exp, condition = pdata[[snakemake@wildcards[["celltype"]]]])
# returns "Default value being used." - this means .5 quantile is being used because the calculated value is < .5
p <- cumNormStatFast(exp)
# fitFeatureModel doesn't work for complicated design matrices
# suggests to use fitZip instead - https://support.bioconductor.org/p/87251/
# res <- fitFeatureModel(exp, design)
print(res)
