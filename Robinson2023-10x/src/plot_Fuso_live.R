library(scater)
library(scran)
library(ggplot2)
library(ggpubr)
library(ggsci)

# counts <- read.table("output/P1/genus_PathSeq_All_reads.tsv", sep="\t", header=TRUE, row.names=1)
counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
# counts <- read.table("output/Pt0/GSM3454529/genus_PathSeq_Bacteria_reads.tsv", sep="\t", header=TRUE, row.names=1)
# output/61_species_PathSeq_metadata.tsv
pdata <- read.table(snakemake@input[[2]], sep="\t", header = TRUE, row.names=1)
row.names(pdata) <- gsub("-", ".", row.names(pdata))
# pdata <- read.table("output/Pt0/GSM3454529/genus_PathSeq_Bacteria_metadata.tsv", sep="\t", header = TRUE, row.names=1)
tax.map <- read.table(snakemake@input[[3]], sep="\t", header=TRUE)
# tax.map <- read.table("output/P1/tax_id_map_All_PathSeq.tsv", sep="\t", header=TRUE)

column.names <- intersect(unique(colnames(counts)), unique(row.names(pdata)))
# remove samples with no spike-ins

counts <- counts[, column.names]
pdata <- pdata[column.names,]

sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData=pdata)

rownames(sce) <- lapply(rownames(sce), function(x) tax.map[tax.map$tax_id == x, "name"])

celltype.col <- "celltype1"
# celltype.col <- snakemake@wildcards[["celltype"]]
# celltype.of.interest <- "HCT116"
celltype.of.interest <- snakemake@wildcards[["celltype_of_interest"]]
# celltype.comparison <- "Jurkat"
celltype.comparison <- snakemake@wildcards[["celltype_comparison"]]
microbe.of.interest <- "Fusobacterium"
min.umis <- strtoi(snakemake@wildcards[["min_umis"]])

# now subset the cells by celltypes of interest if applicable
if (celltype.comparison != "all"){
  sce <- sce[, sce[[celltype.col]] %in% c(celltype.of.interest, celltype.comparison)]
}

celltype <- sce[[celltype.col]]
#
# counts <- t(counts(sce))
# colnames(logcounts) <- lapply(colnames(logcounts), function(x) tax.map[tax.map$tax_id == x, "name"])
# df <- cbind(as.data.frame(t(counts(sce))), as.data.frame(colData(sce)))

get.positive.percent <- function(sce, celltype.of.interest, sample.of.interest) {
  ci.positive <- sum(counts(sce[,(sce[[celltype.col]] == celltype.of.interest) & (sce[["sample"]] == sample.of.interest)])[microbe.of.interest,] >= min.umis)
  ci.total <- dim(sce[,(sce[[celltype.col]] == celltype.of.interest) & (sce[["sample"]] == sample.of.interest)])[[2]]
  ci.percentage <- ci.positive/ci.total
  return(ci.percentage)
}

percent.hct116.3 <- get.positive.percent(sce, "HCT116", "SCAF2963_3_Live")
percent.jurkat.3 <- get.positive.percent(sce, "Jurkat", "SCAF2963_3_Live")
percent.hct116.5 <- get.positive.percent(sce, "HCT116", "SCAF2965_5_Live")
percent.jurkat.5 <- get.positive.percent(sce, "Jurkat", "SCAF2965_5_Live")

# ci.positive <- sum(counts(sce[,(sce[[celltype.col]] == celltype.of.interest) & (sce[["sample"]] == "SCAF2963_3_Live")])[microbe.of.interest,] >= min.umis)
# ci.total <- dim(sce[,celltype == celltype.of.interest])[[2]]
# comparison.positive <- sum(counts(sce[,celltype == celltype.comparison])[microbe.of.interest,] >= min.umis)
# comparison.total <- dim(sce[,celltype == celltype.comparison])[[2]]
# ci.percentage <- ci.positive/ci.total
# comparison.percentage <- comparison.positive/comparison.total
#print(ci.percentage)
#print(comparison.percentage)
df <- data.frame(celltype=c(celltype.comparison, celltype.of.interest, celltype.comparison, celltype.of.interest), percentage=c(percent.hct116.3, percent.jurkat.3, percent.hct116.5, percent.jurkat.5), sample=c("Live_3'", "Live_3'", "Live_5'", "Live_5'"))
df$celltype <- factor(df$celltype, levels=c(celltype.comparison, celltype.of.interest))
#data <- data.frame(group1=c(celltype.of.interest), group2=c(celltype.comparison), p=round(ci.binom[microbe.of.interest, "p.value"], digits=4), y.position=max(ci.percentage, comparison.percentage)*1.1)#

# stat.test <- tibble::tribble(
#   ~group1, ~group2, ~p.adj,
#   "Jurkat", "HCT116", 3.1115424659124e-11,
# )

fig.1d <- ggplot(df, aes(sample, percentage, fill=celltype)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  theme_pubr(base_size=8) + # , base_size= 10 to set font size
  ylab(paste0("% of cells >= ", min.umis, " UMI")) +
  xlab(NULL) +
  # scale_color_npg() +
  scale_fill_manual(values=c("#4DBBD599", "#E64B3599")) + # inspired by scale_color_npg
  scale_y_continuous(labels = scales::percent_format(accuracy = .1)) +
  ggtitle(microbe.of.interest)


# ggbarplot(df, x="celltype", y="percentage")
# ggbarplot(df, aes(celltype, percentage, fill=celltype))
#
# stat_pvalue_manual(
#   stat.test,
#   y.position=.3,
#   label = "p-value"
# )

#
# fig.1d <- ggbarplot(df, x = "celltype", y = "percentage",
#    fill = "celltype") +
#    stat_pvalue_manual(
#      stat.test,
#      y.position=.1
#    ) +
#    scale_fill_manual(values=c("#4DBBD599", "#E64B3599")) +
#    ylab("% of cells >= 2 UMI") +
#    xlab("Celltype") +
#    ggtitle("Fusobacterium") +
#    theme_pubr(legend="none") +
#    scale_y_continuous(labels = scales::percent_format(accuracy = 1))

ggsave(snakemake@output[[1]], width = 2, height = 3, units = "in")
