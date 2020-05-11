library(dplyr)
library(scater)
library(scran)
library(ggplot2)
library(ggpubr)

# use below link as guide
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/

genus.counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
species.counts <- read.table(snakemake@input[[2]], sep="\t", header=TRUE, row.names=1)
# counts <- read.table("output/Pt0_genus_PathSeq_microbe_reads.tsv", sep="\t", header=TRUE, row.names=1)
genus.pdata <- read.table(snakemake@input[[3]], sep="\t", header = TRUE, row.names=1)
species.pdata <- read.table(snakemake@input[[3]], sep="\t", header = TRUE, row.names=1)

# plot genus level finding
min.reads <- 2
min.cells <- 2
genus.counts <- genus.counts[rowSums(genus.counts > min.reads) > min.cells,]
celltype <- genus.pdata[["celltype1"]]
genus.sce <- SingleCellExperiment(assays = list(counts = as.matrix(genus.counts)), colData=genus.pdata)
genus.sce <- computeSumFactors(genus.sce, clusters=celltype)
genus.sce <- logNormCounts(genus.sce)

# Figure 2A1 - expression of Steptomyces
fig.2A1 <- plotExpression(genus.sce, x="celltype1",
    features="Streptomyces") +
    theme(axis.text.x = element_text(angle = 90))

# plot species level findings
min.reads <- 2
min.cells <- 2
species.counts <- species.counts[rowSums(species.counts > min.reads) > min.cells,]
celltype <- species.pdata[["celltype1"]]
species.sce <- SingleCellExperiment(assays = list(counts = as.matrix(species.counts)), colData=species.pdata)
species.sce <- computeSumFactors(species.sce, clusters=celltype)
species.sce <- logNormCounts(species.sce)

# Figure 2A2 - expression of Streptomyces novaecaesareae
fig.2A2 <- plotExpression(species.sce, x="celltype1",
    features="Streptomyces_novaecaesareae") +
    theme(axis.text.x = element_text(angle = 90))

# Figure 2A3 - expression of Streptomyces specialis
fig.2A3 <- plotExpression(species.sce, x="celltype1",
    features="Streptomyces_specialis") +
    theme(axis.text.x = element_text(angle = 90))

ggarrange(fig.2A1, fig.2A2, fig.2A3, ncol=3, nrow=1, labels=c("A", "", ""))

ggsave(snakemake@output[[1]], width=11, height=7)
