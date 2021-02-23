library(scater)
library(scran)
library(ggplot2)
library(ggpubr)
library(ggforce)
# use below link as guide
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/


# read in the read counts files from SRPRISM
counts <- read.table(snakemake@input[[1]], sep="\t", header=FALSE, col.names=c("cell", "reads"), row.names=1)
#counts <- read.table("data/SRPRISM/Pt0/S0/protein_coding_read_counts.tsv", sep="\t", header=FALSE, col.names=c("cell", "reads"), row.names=1)
pdata <- read.table(snakemake@input[[2]], sep="\t", header = TRUE, row.names=2)
#pdata <- read.table("data/units.tsv", sep="\t", header=TRUE, row.names=2)

# subset so both tables have the same rows
row.names <- intersect(unique(rownames(pdata)), unique(rownames(counts)))
pdata <- pdata[row.names,]
counts <- counts[row.names,]
names(counts) <- row.names
# combine them together
fig1a.table <- cbind(pdata, counts)
fig1a.table$logcounts <- log2(fig1a.table$counts + 1)
fig1a.table$infected <- factor(fig1a.table$infected, levels=c("uninfected", "bystander", "infected"))
my.comparisons <- list( c("uninfected", "bystander"), c("bystander", "infected"), c("uninfected", "infected"))


# Figure 1C - plot number of reads that map to protein coding regions

fig.1a <- ggplot(fig1a.table, aes(x=infected, y=logcounts)) +
  geom_violin(scale="width") +
  geom_sina(aes(x = infected,
                y = logcounts,
                fill = infected),
            alpha=0.6,
            size = 2,
            pch=21,
            scale=FALSE) +
  scale_fill_manual(values=c("#FFFF00", "#e74c3c", "#3498db")) +
  ylab("Log2(reads + 1)") +
  xlab("Status") +
  ggtitle("Reads mapped to protein coding regions") +
  theme_pubr(legend="none") +
  stat_compare_means(method="wilcox.test", comparisons = my.comparisons)


ggsave(snakemake@output[[1]])
