library(scater)
library(scran)
library(ggplot2)
library(ggpubr)
library(ggsci)


counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
# counts <- read.table("output/P1/genus_PathSeq_All_reads.tsv", sep="\t", header=TRUE, row.names=1)
# output/61_species_PathSeq_metadata.tsv
pdata <- read.table(snakemake@input[[2]], sep="\t", header = TRUE, row.names=1)
row.names(pdata) <- gsub("-", ".", row.names(pdata))
# pdata <- read.table("output/P1/PathSeq_metadata.tsv", sep="\t", header = TRUE, row.names=1)
tax.map <- read.table(snakemake@input[[3]], sep="\t", header=TRUE)
# tax.map <- read.table("output/P1/tax_id_map_All_PathSeq.tsv", sep="\t", header=TRUE)

column.names <- intersect(unique(colnames(counts)), unique(row.names(pdata)))
# remove samples with no spike-ins

counts <- counts[, column.names]
pdata <- pdata[column.names,]

sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData=pdata)

rownames(sce) <- lapply(rownames(sce), function(x) tax.map[tax.map$tax_id == x, "name"])

microbe.of.interest <- "Fusobacterium"
min.umis <- 1 #strtoi(snakemake@wildcards[["min_umis"]])

get.positive.percent <- function(sce, sample.of.interest) {
  ci.positive <- sum(counts(sce[,(sce[["sample"]] == sample.of.interest)])[microbe.of.interest,] >= min.umis)
  ci.total <- dim(sce[,(sce[["sample"]] == sample.of.interest)])[[2]]
  ci.percentage <- ci.positive/ci.total
  return(ci.percentage)
}


uninfected <- get.positive.percent(sce, "SCAF2961_1_Uninfected")
HK <- get.positive.percent(sce, "SCAF2962_2_HK")
live.3p <- get.positive.percent(sce, "SCAF2963_3_Live")
live.5p <- get.positive.percent(sce, "SCAF2965_5_Live")
condition <- c(rep("Uninfected", 1), rep("Heat-Killed Fn", 1), rep("Live Fn (3' v3)", 1), rep("Live Fn (5')", 1))
positive.percent <- c(uninfected, HK, live.3p, live.5p)

df <- data.frame(condition, positive.percent)

df$condition <- factor(df$condition, levels=c("Uninfected", "Heat-Killed Fn", "Live Fn (3' v3)", "Live Fn (5')"))

ggplot(df, aes(fill=condition, y=positive.percent, x=condition)) +
    geom_bar(position="dodge", stat="identity", show.legend = FALSE) +
    theme_pubr(base_size=6) + # , base_size= 10 to set font size
    # theme(plot.title = element_text(size = 20), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
    ylab("% of cells with >= 1 UMI from Fusobacterium") +
    xlab(NULL) +
    scale_fill_manual(values=c("#00a087", "#4DBBD5FF","#E64B35FF", "#3c5488")) + # inspired by scale_color_npg
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggtitle("in vitro F. nucleatum infection (Robinson2023)")


ggsave(snakemake@output[[1]], width = 3, height = 2, units = "in")
