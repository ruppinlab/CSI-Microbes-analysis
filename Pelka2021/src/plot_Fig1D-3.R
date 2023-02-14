library(scater)
library(scran)
library(ggplot2)
library(ggpubr)
library(ggsci)


# counts <- read.table("output/cohort_reads_for_plot.tsv", sep="\t", header=TRUE, row.names=1)
counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
counts <- t(counts)

# pdata <- read.table("output/cohort_metadata_for_plot.tsv", sep="\t", header = TRUE, row.names=1)
pdata <- read.table(snakemake@input[[2]], sep="\t", header = TRUE, row.names=1)
# row.names(pdata) <- gsub("-", ".", row.names(pdata))


column.names <- intersect(unique(colnames(counts)), unique(row.names(pdata)))
# remove samples with no spike-ins

counts <- counts[, column.names]
pdata <- pdata[column.names,]

sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData=pdata)

# microbe.of.interest <- "Fusobacterium"
min.umis <- 1 #strtoi(snakemake@wildcards[["min_umis"]])

get.positive.percent <- function(sce, chemistry.of.interest, microbe.of.interest) {
  ci.positive <- sum(counts(sce[,(sce[["chemistry"]] == chemistry.of.interest)])[microbe.of.interest,] >= min.umis)
  ci.total <- dim(sce[,(sce[["chemistry"]] == chemistry.of.interest)])[[2]]
  ci.percentage <- ci.positive/ci.total
  return(ci.percentage)
}

get.positive.percent.any <- function(sce, chemistry.of.interest) {
  counts <- counts(sce[,(sce[["chemistry"]] == chemistry.of.interest)])
  ci.positive <- sum(colSums(counts) >= 1)
  ci.total <- dim(sce[,(sce[["chemistry"]] == chemistry.of.interest)])[[2]]
  ci.percentage <- ci.positive/ci.total
  return(ci.percentage)
}


v2.infected.all <- get.positive.percent.any(sce, "v2")
v3.infected.all <- get.positive.percent.any(sce, "v3")
v2.infected.Fusobacterium <- get.positive.percent(sce, "v2", "Fusobacterium")
v3.infected.Fusobacterium <- get.positive.percent(sce, "v3", "Fusobacterium")
v2.infected.Bacteroides <- get.positive.percent(sce, "v2", "Bacteroides")
v3.infected.Bacteroides <- get.positive.percent(sce, "v3", "Bacteroides")


condition <- c(rep("All Genera", 2), rep("Fusobacterium", 2), rep("Bacteroides", 2))
chemistry <- c(rep(c("3' v2", "3' v3"), 3))
percent <- c(v2.infected.all, v3.infected.all, v2.infected.Fusobacterium, v3.infected.Fusobacterium, v2.infected.Bacteroides, v3.infected.Bacteroides)
# percent <- c(1460/1776107498, 9580/1866482021, 41/1776107498, 2732/1866482021, 146/1776107498, 2732/1866482021)

# mat <- matrix(c(9580, 1460, 1866482021, 1776107498), nrow=2, ncol=2, byrow=TRUE)
# f.test <- fisher.test(mat)
# print(f.test)
#
# mat <- matrix(c(2732, 41, 1866482021, 1776107498), nrow=2, ncol=2, byrow=TRUE)
# f.test <- fisher.test(mat)
# print(f.test)
#
# mat <- matrix(c(2732, 146, 1866482021, 1776107498), nrow=2, ncol=2, byrow=TRUE)
# f.test <- fisher.test(mat)
# print(f.test)

df <- data.frame(condition, chemistry, percent)

df$genera <- factor(df$condition, levels=c("All Genera", "Fusobacterium", "Bacteroides"))
df$chemistry <- factor(df$chemistry, levels=c("3' v2", "3' v3"))

ggplot(df, aes(fill=chemistry, y=percent, x=genera)) +
    geom_bar(position="dodge", stat="identity") +
    theme_pubr(base_size=6) + # , base_size= 10 to set font size
    # theme(plot.title = element_text(size = 20), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
    ylab("% of cells with >= 1 UMI") +
    xlab(NULL) +
    scale_fill_manual(values=c("#4DBBD5FF", "#E64B35FF")) + # inspired by scale_color_npg
    scale_y_continuous(labels = scales::percent_format(accuracy = .1)) + # , limits=c(0, .0075)
    ggtitle("Colorectal Carcinoma (Pelka2021)")

# scale_fill_manual() + # inspired by scale_color_npg - values=c("#00A087FF", "#4DBBD5FF","#E64B35FF")

ggsave(snakemake@output[[1]], width = 2.5, height = 2, units = "in")
