library(ggplot2)
library(ggpubr)
library(ggsci)

# use below link as guide
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/


# read in the read counts files from SRPRISM
pdata <- read.table(snakemake@input[[1]], sep="\t", header = TRUE, row.names=2)
#pdata <- read.table("data/units.tsv", sep="\t", header=TRUE, row.names=2)
gene.counts <- read.table(snakemake@input[[2]], sep="\t", header=FALSE, col.names=c("cell", "reads"), row.names=1)
# gene.counts <- read.table("output/Pt0/S0/D23580-gene-read-count.tsv", sep="\t", header=FALSE, col.names=c("cell", "reads"), row.names=1)
rRNA.counts <- read.table(snakemake@input[[3]], sep="\t", header=FALSE, col.names=c("cell", "reads"), row.names=1)
# rRNA.counts <- read.table("output/Pt0/S0/D23580-rRNA-read-count.tsv", sep="\t", header=FALSE, col.names=c("cell", "reads"), row.names=1)
# print(gene.counts)
# print(rRNA.counts)
counts <- gene.counts + rRNA.counts
# print(counts)
# subset so both tables have the same rows
row.names <- intersect(unique(rownames(pdata)), unique(rownames(counts)))
pdata <- pdata[row.names,]
counts <- counts[row.names,]
names(counts) <- row.names

# combine them together
df <- cbind(pdata, counts)
df$logcounts <- log2(df$counts + 1)
df$infected <- factor(df$infected, levels=c("uninfected", "bystander", "infected"))
levels(df$infected) <- c("Control", "Bystander", "Infected")

#print(df[df$infected == "Infected", "logcounts"])
print(wilcox.test(df[df$infected == "Infected", "logcounts"], df[df$infected == "Bystander", "logcounts"]))
print(wilcox.test(df[df$infected == "Infected", "logcounts"], df[df$infected == "Control", "logcounts"]))
# Figure 1A - plot total number of reads
#my.comparisons <- list( c("uninfected", "bystander"), c("bystander", "infected"), c("uninfected", "infected"))
# ylabel <- bquote("" * log[2] * "(reads + 1)")
#ylabel <- expression(log[2](reads+1))
#print(ylabel)
# fig.1a <- ggboxplot(fig1a.table, x="infected", y="logcounts", color="infected",
#   add="jitter", palette = c("#FFFF00", "#e74c3c", "#3498db")
# )

# geom_jitter() is shortcut for geom_point(position = "jitter")
fig.1a <- ggplot(df, aes(x=infected, y=logcounts)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(aes(color=infected), alpha=0.6, size=2) +
  theme_pubr(legend="none", base_size=9) + # , base_size= 10 to set font size
  ylab("log2(reads+1)") +
  xlab(NULL) +
  ggtitle("Reads mapped to \nSalmonella genome") +
  #scale_color_npg()
  scale_color_manual(values=c("#00A087FF", "#4DBBD5FF","#E64B35FF")) # inspired by scale_color_npg
  #theme(plot.title = element_text(size=12), axis.title.y = element_text(size=10))
  #stat_compare_means(method="wilcox.test", comparisons = my.comparisons, size=1) +

  #
  # geom_sina(aes(x = infected,
  #                 y = logcounts,
  #                 fill = infected),
  #             alpha=0.6,
  #             size = 2,
  #             pch=21,
  #             scale=FALSE) +

# fig.1a <- ggplot(fig1a.table, aes(x=infected, y=logcounts)) +
#   geom_violin(scale="width") +
#   geom_sina(aes(x = infected,
#                 y = logcounts,
#                 fill = infected),
#             alpha=0.6,
#             size = 2,
#             pch=21,
#             scale=FALSE) +
#   scale_fill_manual(values=c("#FFFF00", "#e74c3c", "#3498db")) +
#   ylab("Log2(reads + 1)") +
#   # xlab("Status") +
#   ggtitle("Reads mapped to entire Salmonella genome") +
#   theme_pubr(legend="none") +
#   stat_compare_means(method="wilcox.test", comparisons = my.comparisons)

# need to change the size here
ggsave(snakemake@output[[1]], width = 3.25, height = 3.5, units = "in")
