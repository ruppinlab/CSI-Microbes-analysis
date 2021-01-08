library(scater)
library(scran)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(stringr)
# use below link as guide
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/

source(snakemake@params[["spike_functions"]])

sce <- generate_sce(microbe.file=snakemake@input[[1]], spikein.file=snakemake@input[[2]],
  pdata.file=snakemake@input[[3]], celltype.col="celltype1")

tax.map <- read.table(snakemake@input[[4]], sep="\t", header=TRUE)
tax.map <- tax.map[tax.map$taxa_level == "species",]

# rename rownames from tax_id to names
rownames(sce) <- lapply(rownames(sce), function(x) tax.map[tax.map$tax_id == x, "name"])

sce <- sce[, sce[["celltype1"]] %in% c("Tumor", "epithelial", "immune")]

plot_fig1 <- data.frame(
    logcounts = logcounts(sce)["Cutibacterium_acnes", ],
    status = sce[["celltype1"]])
print(plot_fig1)
fig.1a <- ggplot(plot_fig1, aes(x=status, y=logcounts)) +
geom_violin(scale="width") +
geom_sina(aes(x = status,
  y = logcounts,
  fill = status),
  alpha=0.6,
  size = 2,
  pch=21,
  scale=FALSE) +
  scale_fill_manual(values=c("#e74c3c", "#3498db", "008000")) +
  ylab("Expression log(counts)") +
  xlab("Cells") +
  ggtitle("Cutibacterium acnes") +
  theme_pubr(legend="none")

sce <- generate_sce(microbe.file=snakemake@input[[5]], spikein.file=snakemake@input[[6]],
  pdata.file=snakemake@input[[7]], celltype.col="celltype1")

tax.map <- read.table(snakemake@input[[8]], sep="\t", header=TRUE)
tax.map <- tax.map[tax.map$taxa_level == "species",]

# rename rownames from tax_id to names
rownames(sce) <- lapply(rownames(sce), function(x) tax.map[tax.map$tax_id == x, "name"])

sce <- sce[, sce[["celltype1"]] %in% c("Tumor", "immune", "epithelial")]

plot_fig1 <- data.frame(
    logcounts = logcounts(sce)["Gardnerella_vaginalis", ],
    status = sce[["celltype1"]])
print(plot_fig1)
fig.1b <- ggplot(plot_fig1, aes(x=status, y=logcounts)) +
geom_violin(scale="width") +
geom_sina(aes(x = status,
  y = logcounts,
  fill = status),
  alpha=0.6,
  size = 2,
  pch=21,
  scale=FALSE) +
  scale_fill_manual(values=c("#e74c3c", "#3498db", "008000")) +
  ylab("Expression log(counts)") +
  xlab("Cells") +
  ggtitle("Gardnerella vaginalis") +
  theme_pubr(legend="none")



ggarrange(fig.1a, fig.1b, ncol=2, nrow=1, labels=c("A", "B"))

ggsave(snakemake@output[[1]], width=16, height=6)
