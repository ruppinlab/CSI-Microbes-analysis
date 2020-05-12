library(dplyr)
library(scater)
library(scran)
library(ggplot2)
library(ggpubr)

# use below link as guide
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/

mel129pa.counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
mel129pa.pdata <- read.table(snakemake@input[[2]], sep="\t", header=TRUE, row.names=1)
# counts <- read.table("output/Pt0_genus_PathSeq_microbe_reads.tsv", sep="\t", header=TRUE, row.names=1)
mel106.counts <- read.table(snakemake@input[[3]], sep="\t", header = TRUE, row.names=1)
mel106.pdata <- read.table(snakemake@input[[4]], sep="\t", header = TRUE, row.names=1)

# plot mel129pa findings
min.reads <- 2
min.cells <- 2
mel129pa.counts <- mel129pa.counts[rowSums(mel129pa.counts > min.reads) > min.cells,]
celltype <- mel129pa.pdata[["Tumor"]]
mel129pa.sce <- SingleCellExperiment(assays = list(counts = as.matrix(mel129pa.counts)), colData=mel129pa.pdata)
mel129pa.sce <- computeSumFactors(mel129pa.sce, clusters=celltype)
mel129pa.sce <- logNormCounts(mel129pa.sce)

# Figure 2C1 - expression of Propionibacterium
fig.2C1 <- plotExpression(mel129pa.sce, x="celltype1",
    features="Propionibacterium") +
    theme(axis.text.x = element_text(angle = 90))

# plot species level findings
min.reads <- 2
min.cells <- 2
mel106.counts <- mel106.counts[rowSums(mel106.counts > min.reads) > min.cells,]
celltype <- mel106.pdata[["celltype1"]]
mel106.sce <- SingleCellExperiment(assays = list(counts = as.matrix(mel106.counts)), colData=mel106.pdata)
mel106.sce <- computeSumFactors(mel106.sce, clusters=celltype)
mel106.sce <- logNormCounts(mel106.sce)

# Figure 2C2 - expression of Coprinopsis_cinerea
fig.2C2 <- plotExpression(mel106.sce, x="celltype1",
    features="Coprinopsis_cinerea") +
    theme(axis.text.x = element_text(angle = 90))


ggarrange(fig.2C1, fig.2C2, ncol=2, nrow=1, labels=c("C", "", ""))

ggsave(snakemake@output[[1]], width=11, height=7)
