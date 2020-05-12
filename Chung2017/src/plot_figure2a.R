library(dplyr)
library(scater)
library(scran)
library(ggplot2)
library(ggforce)
library(ggpubr)

# use below link as guide
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/

setwd("~/src/CSI-Microbes-analysis/Chung2017/")

#genus.counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
#species.counts <- read.table(snakemake@input[[2]], sep="\t", header=TRUE, row.names=1)
genus.counts <- read.table("output/BC06_genus_PathSeq_microbe_reads.tsv", sep="\t", header=TRUE, row.names=1)
species.counts <- read.table("output/BC06_species_PathSeq_microbe_reads.tsv", sep="\t", header=TRUE, row.names=1)
# counts <- read.table("output/Pt0_genus_PathSeq_microbe_reads.tsv", sep="\t", header=TRUE, row.names=1)
#genus.pdata <- read.table(snakemake@input[[3]], sep="\t", header = TRUE, row.names=1)
#species.pdata <- read.table(snakemake@input[[3]], sep="\t", header = TRUE, row.names=1)
genus.pdata <- read.table("output/BC06_genus_PathSeq_metadata.tsv", sep="\t", header=TRUE, row.names=1)
species.pdata <- read.table("output/BC06_species_PathSeq_metadata.tsv", sep="\t", header=TRUE, row.names=1)

# plot genus level finding
min.reads <- 2
min.cells <- 2
genus.counts <- genus.counts[rowSums(genus.counts > min.reads) > min.cells,]
celltype <- genus.pdata[["celltype1"]]
genus.sce <- SingleCellExperiment(assays = list(counts = as.matrix(genus.counts)), colData=genus.pdata)
genus.sce <- computeSumFactors(genus.sce, clusters=celltype)
genus.sce <- logNormCounts(genus.sce)

# Figure 2A1 - expression of Steptomyces
feature <- "Streptomyces"
sce_logcounts <- genus.sce@assays@data$logcounts
plot_fig2 <- data.frame(
    logcounts = sce_logcounts[feature, ],
    cell_types = genus.sce$celltype1)

fig.2A1 <- ggplot(plot_fig2, aes(x=cell_types, y=logcounts)) + 
    geom_violin(scale="width", adjust=0.7, trim=TRUE) +
    geom_sina(aes(x = cell_types, 
                  y = logcounts, 
                  fill = cell_types),
              alpha=0.6,
              size = 3,
              pch=21,
              scale=FALSE) +
    ylim(c(0,13)) +
    scale_fill_manual(values=c("#e74c3c", "#3498db", "#2ecc71")) +
    ylab("Expression log(counts)") +
    xlab("Cells") +
    ggtitle(feature) +
    theme_pubr(legend="none")
plot(fig.2A1)

# plot species level findings
min.reads <- 2
min.cells <- 2
species.counts <- species.counts[rowSums(species.counts > min.reads) > min.cells,]
celltype <- species.pdata[["celltype1"]]
species.sce <- SingleCellExperiment(assays = list(counts = as.matrix(species.counts)), colData=species.pdata)
species.sce <- computeSumFactors(species.sce, clusters=celltype)
species.sce <- logNormCounts(species.sce)

# Figure 2A2 - expression of Streptomyces novaecaesareae
feature <- "Streptomyces_novaecaesareae"
sce_logcounts <- species.sce@assays@data$logcounts
plot_fig2 <- data.frame(
    logcounts = sce_logcounts[feature, ],
    cell_types = genus.sce$celltype1)

fig.2A2 <- ggplot(plot_fig2, aes(x=cell_types, y=logcounts)) + 
    geom_violin(scale="width", trim=TRUE) +
    geom_jitter(aes(x=cell_types, y=logcounts, fill=celltype), pch=21, alpha=0.6, size=3) +
    scale_fill_manual(values=c("#e74c3c", "#3498db", "#2ecc71")) +
    ylab("Expression log(counts)") +
    xlab("Cells") +
    ggtitle("Streptomyces novaecaesareae") +
    theme_pubr(legend="none")
plot(fig.2A2)

# Figure 2A3 - expression of Streptomyces specialis
feature <- "Streptomyces_specialis"
sce_logcounts <- species.sce@assays@data$logcounts
plot_fig2 <- data.frame(
    logcounts = sce_logcounts[feature, ],
    cell_types = genus.sce$celltype1)

fig.2A3 <- ggplot(plot_fig2, aes(x=cell_types, y=logcounts)) + 
    geom_violin(scale="width", trim=TRUE) +
    geom_jitter(aes(x=cell_types, y=logcounts, fill=celltype), pch=21, alpha=0.6, size=3) +
    scale_fill_manual(values=c("#e74c3c", "#3498db", "#2ecc71")) +
    ylab("Expression log(counts)") +
    xlab("Cells") +
    ggtitle("Streptomyces specialis") +
    theme_pubr(legend="none")
plot(fig.2A3)

ggarrange(fig.2A1, fig.2A2, fig.2A3, ncol=3, nrow=1, labels=c("A", "", ""))
#ggsave(filename="fig2.pdf", path="~/src/CSI-Microbes-analysis/Chung2017/output/", width=11, height=5)

ggsave(snakemake@output[[1]], width=11, height=5)
