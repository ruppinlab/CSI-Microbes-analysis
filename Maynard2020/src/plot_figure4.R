library(scater)
library(scran)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(stringr)
# use below link as guide
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/

source(snakemake@params[["spike_functions"]])

# Plot Cutibacterium acnes for TH236
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
  scale_fill_manual(values=c("#e74c3c", "#3498db", "#000000")) +
  ylab("Expression log(counts)") +
  xlab("Cells") +
  ggtitle("Cutibacterium acnes abundance in patient TH236") +
  theme_pubr(legend="none")


# plot Corynebacterium for TH238
sce <- generate_sce(microbe.file=snakemake@input[[5]], spikein.file=snakemake@input[[6]],
  pdata.file=snakemake@input[[7]], celltype.col="celltype1")
print(sce)
tax.map <- read.table(snakemake@input[[8]], sep="\t", header=TRUE)
tax.map <- tax.map[tax.map$taxa_level == "genus",]

# rename rownames from tax_id to names
rownames(sce) <- lapply(rownames(sce), function(x) tax.map[tax.map$tax_id == x, "name"])

#sce <- sce[, sce[["celltype1"]] %in% c("Tumor", "immune", "epithelial")]

plot_fig1 <- data.frame(
    logcounts = logcounts(sce)["Corynebacterium", ],
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
  scale_fill_manual(values=c("#e74c3c", "#3498db", "008000", "#000000")) +
  ylab("Expression log(counts)") +
  xlab("Cells") +
  ggtitle("Corynebacterium abundance in patient TH238") +
  theme_pubr(legend="none")

# plot Cutibacterium acnes for TH266
sce <- generate_sce(microbe.file=snakemake@input[[9]], spikein.file=snakemake@input[[10]],
  pdata.file=snakemake@input[[11]], celltype.col="celltype1")

tax.map <- read.table(snakemake@input[[12]], sep="\t", header=TRUE)
tax.map <- tax.map[tax.map$taxa_level == "species",]

# rename rownames from tax_id to names
rownames(sce) <- lapply(rownames(sce), function(x) tax.map[tax.map$tax_id == x, "name"])

#sce <- sce[, sce[["celltype1"]] %in% c("Tumor", "immune", "stroma")]

plot_fig1 <- data.frame(
    logcounts = logcounts(sce)["Cutibacterium_acnes", ],
    status = sce[["celltype1"]])
print(plot_fig1)
fig.1c <- ggplot(plot_fig1, aes(x=status, y=logcounts)) +
geom_violin(scale="width") +
geom_sina(aes(x = status,
  y = logcounts,
  fill = status),
  alpha=0.6,
  size = 2,
  pch=21,
  scale=FALSE) +
  scale_fill_manual(values=c("#3498db", "008000", "#000000")) +
  ylab("Expression log(counts)") +
  xlab("Cells") +
  ggtitle("Cutibacterium acnes abundance in patient TH266") +
  theme_pubr(legend="none")

  # plot Leptotrichia for TH231
  sce <- generate_sce(microbe.file=snakemake@input[[13]], spikein.file=snakemake@input[[14]],
    pdata.file=snakemake@input[[15]], celltype.col="celltype1")

  tax.map <- read.table(snakemake@input[[16]], sep="\t", header=TRUE)
  tax.map <- tax.map[tax.map$taxa_level == "genus",]

  # rename rownames from tax_id to names
  rownames(sce) <- lapply(rownames(sce), function(x) tax.map[tax.map$tax_id == x, "name"])

  sce <- sce[, sce[["celltype1"]] %in% c("Tumor", "immune", "stroma")]

  plot_fig1 <- data.frame(
      logcounts = logcounts(sce)["Leptotrichia", ],
      status = sce[["celltype1"]])
  print(plot_fig1)
  fig.1d <- ggplot(plot_fig1, aes(x=status, y=logcounts)) +
  geom_violin(scale="width") +
  geom_sina(aes(x = status,
    y = logcounts,
    fill = status),
    alpha=0.6,
    size = 2,
    pch=21,
    scale=FALSE) +
    scale_fill_manual(values=c("#3498db", "008000", "#000000")) +
    ylab("Expression log(counts)") +
    xlab("Cells") +
    ggtitle("Leptotrichia abundance in patient TH231") +
    theme_pubr(legend="none")

ggarrange(fig.1a, fig.1b, fig.1c, fig.1d, ncol=2, nrow=2, labels=c("A", "B", "C", "D"))

ggsave(snakemake@output[[1]], width=16, height=10)
