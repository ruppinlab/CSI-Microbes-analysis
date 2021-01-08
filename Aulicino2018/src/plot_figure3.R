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
  pdata.file=snakemake@input[[3]], celltype.col="infected")

tax.map <- read.table(snakemake@input[[4]], sep="\t", header=TRUE)
tax.map <- tax.map[tax.map$taxa_level == "species",]

# rename rownames from tax_id to names
rownames(sce) <- lapply(rownames(sce), function(x) tax.map[tax.map$tax_id == x, "name"])

sce <- sce[, sce[["infected"]] %in% c("infected", "bystander")]

plot_fig1 <- data.frame(
    logcounts = logcounts(sce)["Salmonella_enterica", ],
    status = sce$infected)
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
  scale_fill_manual(values=c("#e74c3c", "#3498db")) +
  ylab("Expression log(counts)") +
  xlab("Cells") +
  ggtitle("Salmonella enterica") +
  theme_pubr(legend="none")

# stats
block <- sce[["plate"]]
wilcox <- findMarkers(sce, test="wilcox", groups=sce[["infected"]], lfc=0.5, block=block)
df <- wilcox[["infected"]]
df$Species <- row.names(df)
print(df$Species)
df$Species <- gsub("_", " ", df$Species)
df$Species = str_wrap(df$Species, width = 15)
df <- df[0:10,]

df <- as.data.frame(df)
df <- df[order(df$AUC.bystander, decreasing=FALSE), ]
df$Species <- factor(df$Species, levels=df$Species)

# Try to plot
fig.1d <- ggplot(data=df, aes(x=Species, y=AUC.bystander)) +
  geom_bar(stat="identity") +
  geom_text(data=df, aes(
    x=Species,
    y=AUC.bystander-0.2,
    label=paste("FDR pval=", round(FDR, digits = 3) , sep="")),
    colour="white") +
  coord_flip() +
  theme_pubr(legend="bottom") +
  theme(axis.text = element_text(size = 10))


# Plot Class Results
sce <- generate_sce(microbe.file=snakemake@input[[5]], spikein.file=snakemake@input[[2]],
  pdata.file=snakemake@input[[3]], celltype.col="infected")

tax.map <- read.table(snakemake@input[[4]], sep="\t", header=TRUE)
tax.map <- tax.map[tax.map$taxa_level == "class",]

# rename rownames from tax_id to names
rownames(sce) <- lapply(rownames(sce), function(x) tax.map[tax.map$tax_id == x, "name"])

sce <- sce[, sce[["infected"]] %in% c("infected", "bystander")]

plot_fig1 <- data.frame(
    logcounts = logcounts(sce)["Gammaproteobacteria", ],
    status = sce$infected)
#print(plot_fig1)
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
  ggtitle("Gammaproteobacteria") +
  theme_pubr(legend="none")

# stats
block <- sce[["plate"]]
wilcox <- findMarkers(sce, test="wilcox", groups=sce[["infected"]], lfc=0.5, block=block)
df <- wilcox[["infected"]]
df$Classes <- row.names(df)
df <- df[0:10,]

df <- as.data.frame(df)
df <- df[order(df$AUC.bystander, decreasing=FALSE), ]
df$Classes <- factor(df$Classes, levels=df$Classes)

# Try to plot
fig.1b <- ggplot(data=df, aes(x=Classes, y=AUC.bystander)) +
  geom_bar(stat="identity") +
  geom_text(data=df, aes(
    x=Classes,
    y=AUC.bystander-0.2,
    label=paste("FDR pval=", round(FDR, digits = 3) , sep="")),
    colour="white") +
  coord_flip() +
  theme_pubr(legend="bottom") +
  theme(axis.text = element_text(size = 10))
ggarrange(fig.1a, fig.1b, fig.1c, fig.1d, ncol=4, nrow=1, labels=c("A", "B", "C", "D"))

ggsave(snakemake@output[[1]], width=16, height=6)
