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

get.positive.percent <- function(sce, celltype.of.interest, sample.of.interest) {
  ci.positive <- sum(counts(sce[,(sce[["celltype1"]] == celltype.of.interest) & (sce[["sample"]] == sample.of.interest)])[microbe.of.interest,] >= min.umis)
  ci.total <- dim(sce[,(sce[["celltype1"]] == celltype.of.interest) & (sce[["sample"]] == sample.of.interest)])[[2]]
  ci.percentage <- ci.positive/ci.total
  return(ci.percentage)
}


Jurkat.uninfected <- get.positive.percent(sce, "Jurkat", "SCAF2961_1_Uninfected")
HCT116.uninfected <- get.positive.percent(sce, "HCT116", "SCAF2961_1_Uninfected")
THP1.uninfected <- get.positive.percent(sce, "THP1", "SCAF2961_1_Uninfected")
Jurkat.HK <- get.positive.percent(sce, "Jurkat", "SCAF2962_2_HK")
HCT116.HK <- get.positive.percent(sce, "HCT116", "SCAF2962_2_HK")
THP1.HK <- get.positive.percent(sce, "THP1", "SCAF2962_2_HK")
Jurkat.live.3p <- get.positive.percent(sce, "Jurkat", "SCAF2963_3_Live")
HCT116.live.3p <- get.positive.percent(sce, "HCT116", "SCAF2963_3_Live")
THP1.live.3p <- get.positive.percent(sce, "THP1", "SCAF2963_3_Live")
Jurkat.live.5p <- get.positive.percent(sce, "Jurkat", "SCAF2965_5_Live")
HCT116.live.5p <- get.positive.percent(sce, "HCT116", "SCAF2965_5_Live")
THP1.live.5p <- get.positive.percent(sce, "THP1", "SCAF2965_5_Live")
condition <- c(rep("Uninfected (3' v3)", 3), rep("Heat-Killed Fn (3' v3)", 3), rep("Live Fn (3' v3)", 3), rep("Live Fn (5')", 3))
celltype <- c(rep(c("Jurkat", "HCT116", "THP1"), 4))
positive.percent <- c(Jurkat.uninfected, HCT116.uninfected, THP1.uninfected, Jurkat.HK, HCT116.HK, THP1.HK, Jurkat.live.3p, HCT116.live.3p, THP1.live.3p, Jurkat.live.5p, HCT116.live.5p, THP1.live.5p)

df <- data.frame(condition, celltype, positive.percent)

df$condition <- factor(df$condition, levels=c("Uninfected (3' v3)", "Heat-Killed Fn (3' v3)", "Live Fn (3' v3)", "Live Fn (5')"))
df$celltype <- factor(df$celltype, levels=c("Jurkat", "HCT116", "THP1"))

ggplot(df, aes(fill=celltype, y=positive.percent, x=condition)) +
    geom_bar(position="dodge", stat="identity") +
    theme_pubr(base_size=6) + # , base_size= 10 to set font size
    # theme(plot.title = element_text(size = 20), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
    ylab("% of cells with >= 1 UMI from Fusobacterium") +
    xlab(NULL) +
    scale_fill_manual(values=c("#4DBBD5FF","#E64B35FF", "#3c5488")) + # inspired by scale_color_npg
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggtitle("in vitro F. nucleatum infection (Robinson2023)")


ggsave(snakemake@output[[1]], width = 3.5, height = 2, units = "in")
