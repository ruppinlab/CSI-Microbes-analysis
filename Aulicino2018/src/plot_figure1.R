library(scater)
library(scran)
library(ggplot2)
library(ggpubr)
library(ggforce)
# use below link as guide
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/

#setwd("~/src/CSI-Microbes-analysis/Aulicino2018/")

#counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
counts <- read.table("output/Pt0_genus_PathSeq_microbe_reads.tsv", sep="\t", header=TRUE, row.names=1)
#pdata <- read.table(snakemake@input[[2]], sep="\t", header = TRUE, row.names=1)
pdata <- read.table("output/Pt0_genus_PathSeq_metadata.tsv", sep="\t", header=TRUE, row.names=1)

row.names(pdata) <- gsub("-", ".", row.names(pdata))
# pdata <- read.table("output/Pt0_PathSeq_metadata.tsv", sep="\t", header = TRUE, row.names=1)

# filter out any celltypes with < 5 cells
celltype.col <- "infected"
agg <- aggregate(row.names(pdata), by=list(pdata[[celltype.col]]), FUN=length)
celltypes.to.keep <- agg[agg$x >= 5,]$Group.1
pdata <- pdata[pdata[[celltype.col]] %in% celltypes.to.keep, ]
# filter out unknown celltypes
pdata <- pdata[pdata[[celltype.col]] != "unknown", ]
pdata <- droplevels(pdata)
counts <- counts[, row.names(pdata)]

# filter out any OTUs without the minimum number of reads in .7 of the smallest group
# before normalization -> reduce sparsity to manageable amount
min.reads <- 2
min.cells <- min(agg$x)*.7
counts <- counts[rowSums(counts > min.reads) > min.cells,]
# pdata <- pdata[pdata["selection"] == "Astrocytes(HEPACAM)", ]
# pdata <- pdata[pdata["depletion_batch"] != "depleted_yes", ]
celltype <- pdata[[celltype.col]]
sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData=pdata)
sce <- computeSumFactors(sce, clusters=celltype)
sce <- logNormCounts(sce)

# Figure 1A - plot expression
# Extract logcounts from single cell experiment
feature <- "Salmonella"
sce_logcounts <- sce@assays@data$logcounts
plot_fig1 <- data.frame(
    logcounts = sce_logcounts[feature, ],
    status = sce$infected)

fig.1a <- ggplot(plot_fig1, aes(x=status, y=logcounts)) + 
  geom_violin(scale="width") +
  geom_sina(aes(x = status, 
                y = logcounts, 
                fill = status),
            alpha=0.6,
            size = 2,
            pch=21,
            scale=FALSE) +
  scale_fill_manual(values=c("#e74c3c", "#3498db")) +
  ylab("Expression log(counts)") +
  xlab("Cells") +
  ggtitle(feature) +
  theme_pubr(legend="bottom")

# stats
wilcox <- findMarkers(sce, test="wilcox", groups=celltype, lfc=1, block=pdata$plate)
df <- wilcox[["infected"]]
df$Genera <- row.names(df)
df <- df[0:25,]
print(df)

# Try to plot

# Figure 1B - plot table
fig.1b <- ggtexttable(df[, c("Genera", "FDR", "AUC.uninfected")], rows=NULL)

ggarrange(fig.1a, fig.1b, ncol=2, nrow=1, labels=c("B", "C"))

ggsave(snakemake@output[[1]], width=11, height=7)
