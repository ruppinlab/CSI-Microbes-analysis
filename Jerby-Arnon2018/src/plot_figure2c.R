library(dplyr)
library(scater)
library(scran)
library(ggplot2)
library(ggpubr)
library(ggforce)

# use below link as guide
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/

#setwd("~/src/CSI-Microbes-analysis/Jerby-Arnon2018/")

mel129pa.counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
mel129pa.pdata <- read.table(snakemake@input[[2]], sep="\t", header=TRUE, row.names=1)
# counts <- read.table("output/Pt0_genus_PathSeq_microbe_reads.tsv", sep="\t", header=TRUE, row.names=1)
mel106.counts <- read.table(snakemake@input[[3]], sep="\t", header = TRUE, row.names=1)
mel106.pdata <- read.table(snakemake@input[[4]], sep="\t", header = TRUE, row.names=1)

#mel129pa.counts <- read.table("output/Mel129pa_genus_PathSeq_microbe_reads.tsv", sep="\t", header=TRUE, row.names=1)
#mel129pa.pdata <- read.table("output/Mel129pa_genus_PathSeq_metadata.tsv", sep="\t", header=TRUE, row.names=1)
#mel106.counts <- read.table("output/Mel106_species_PathSeq_microbe_reads.tsv", sep="\t", header=TRUE, row.names=1)
#mel106.pdata <- read.table("output/Mel106_species_PathSeq_metadata.tsv", sep="\t", header=TRUE, row.names=1)

# plot mel129pa findings
min.reads <- 2
min.cells <- 2
mel129pa.counts <- mel129pa.counts[rowSums(mel129pa.counts > min.reads) > min.cells,]
celltype <- mel129pa.pdata[["Tumor"]]
mel129pa.sce <- SingleCellExperiment(assays = list(counts = as.matrix(mel129pa.counts)), colData=mel129pa.pdata)
mel129pa.sce <- computeSumFactors(mel129pa.sce, clusters=celltype)
mel129pa.sce <- logNormCounts(mel129pa.sce)

# Figure 2B1 - expression of Propionibacterium
feature <- "Propionibacterium"
sce_logcounts <- mel129pa.sce@assays@data$logcounts
plot_fig1 <- data.frame(
  logcounts = sce_logcounts[feature, ],
  status = mel129pa.sce$celltype1)

fig.2B1 <- ggplot(plot_fig1, aes(x=status, y=logcounts)) +
  geom_boxplot(data=plot_fig1, outlier.shape = NA) +
  #geom_violin(scale="width") +
  geom_sina(aes(x = status,
                y = logcounts,
                fill = status),
            alpha=0.6,
            size = 2,
            pch=21,
            scale=FALSE) +
  #scale_fill_manual(values=c("#e74c3c", "#3498db")) +
  ylab("Expression log(counts)") +
  xlab("Cells") +
  ggtitle(feature) +
  theme_pubr(legend="none")

plot(fig.2B1)


# plot species level findings
min.reads <- 2
min.cells <- 2
mel106.counts <- mel106.counts[rowSums(mel106.counts > min.reads) > min.cells,]
celltype <- mel106.pdata[["celltype1"]]
mel106.sce <- SingleCellExperiment(assays = list(counts = as.matrix(mel106.counts)), colData=mel106.pdata)
mel106.sce <- computeSumFactors(mel106.sce, clusters=celltype)
mel106.sce <- logNormCounts(mel106.sce)

# Figure 2B2 - expression of Coprinopsis_cinerea
feature <- "Coprinopsis_cinerea"
sce_logcounts <- mel106.sce@assays@data$logcounts
plot_fig1 <- data.frame(
  logcounts = sce_logcounts[feature, ],
  status = mel106.sce$celltype1)

fig.2B2 <- ggplot(plot_fig1, aes(x=status, y=logcounts)) +
  geom_boxplot(data=plot_fig1, outlier.shape = NA) +
  #geom_violin(scale="width") +
  geom_sina(aes(x = status,
                y = logcounts,
                fill = status),
            alpha=0.6,
            size = 2,
            pch=21,
            scale=FALSE) +
  #scale_fill_manual(values=c("#e74c3c", "#3498db")) +
  ylab("Expression log(counts)") +
  xlab("Cells") +
  ggtitle(feature) +
  theme_pubr(legend="none")

plot(fig.2B2)


ggarrange(fig.2B1, fig.2B2, ncol=2, nrow=1, labels=c("", "", ""))
#ggsave(filename="fig2b_boxplot.pdf", path="~/src/CSI-Microbes-analysis/Jerby-Arnon2018/", width=9, height=4)


ggsave(snakemake@output[[1]], width=11, height=7)
