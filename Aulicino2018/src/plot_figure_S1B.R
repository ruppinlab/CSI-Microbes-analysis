library(ggplot2)
library(ggpubr)
library(ggsci)


# use below link as guide
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/


# read in the read counts files from SRPRISM
counts <- read.table(snakemake@input[[1]], sep="\t", header=FALSE, col.names=c("cell", "reads"), row.names=1)
#counts <- read.table("data/SRPRISM/Pt0/S0/rRNA_read_counts.tsv", sep="\t", header=FALSE, col.names=c("cell", "reads"), row.names=1)
pdata <- read.table(snakemake@input[[2]], sep="\t", header = TRUE, row.names=2)
#pdata <- read.table("data/units.tsv", sep="\t", header=TRUE, row.names=2)
print(counts)
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
#my.comparisons <- list( c("uninfected", "bystander"), c("bystander", "infected"), c("uninfected", "infected"))

print(wilcox.test(df[df$infected == "Infected", "logcounts"], df[df$infected == "Bystander", "logcounts"]))
print(wilcox.test(df[df$infected == "Infected", "logcounts"], df[df$infected == "Control", "logcounts"]))

# Figure 1B - plot number of reads mapped to rRNA operands

fig.1b <- ggplot(df, aes(x=infected, y=logcounts)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(aes(color=infected), alpha=0.6, size=2) +
  theme_pubr(legend="none", base_size=9) + # , base_size= 10 to set font size
  ylab("log2(reads+1)") +
  xlab(NULL) +
  ggtitle("Reads mapped to \nrRNA regions") +
  #scale_color_npg()
  scale_color_manual(values=c("#00A087FF", "#4DBBD5FF","#E64B35FF")) # inspired by scale_color_npg

ggsave(snakemake@output[[1]], width = 3.25, height = 3.5, units = "in")
