library(scater)
library(scran)
library(ggplot2)
library(ggpubr)
library(ggsci)


counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
# counts <- read.table("output/Pt0/genus_PathSeq_All_reads.tsv", sep="\t", header=TRUE, row.names=1)
# output/61_species_PathSeq_metadata.tsv
pdata <- read.table(snakemake@input[[2]], sep="\t", header = TRUE, row.names=1)
# pdata <- read.table("output/Pt0/PathSeq_metadata.tsv", sep="\t", header = TRUE, row.names=1)
row.names(pdata) <- gsub("-", ".", row.names(pdata))
tax.map <- read.table(snakemake@input[[3]], sep="\t", header=TRUE)
# tax.map <- read.table("output/Pt0/tax_id_map_All_PathSeq.tsv", sep="\t", header=TRUE)

column.names <- intersect(unique(colnames(counts)), unique(row.names(pdata)))
# remove samples with no spike-ins

counts <- counts[, column.names]
pdata <- pdata[column.names,]

sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData=pdata)

rownames(sce) <- lapply(rownames(sce), function(x) tax.map[tax.map$tax_id == x, "name"])

celltype.col <- "Monocyte"
# celltype.col <- snakemake@wildcards[["celltype"]]
#celltype.of.interest <- "Monocyte"
# celltype.of.interest <- snakemake@wildcards[["celltype_of_interest"]]
#celltype.comparison <- "nonMonocyte"
# celltype.comparison <- snakemake@wildcards[["celltype_comparison"]]
#microbe.of.interest <- "Salmonella"
#microbe.of.interest <- snakemake@wildcards[["microbe_of_interest"]]

# now subset the cells by celltypes of interest if applicable
# if (celltype.comparison != "all"){
#   sce <- sce[, sce[[celltype.col]] %in% c(celltype.of.interest, celltype.comparison)]
# }

# celltype <- sce[[celltype.col]]
#
# counts <- t(counts(sce))
# colnames(logcounts) <- lapply(colnames(logcounts), function(x) tax.map[tax.map$tax_id == x, "name"])
# df <- cbind(as.data.frame(t(counts(sce))), as.data.frame(colData(sce)))

min.umis <- 1

get.positive.percent <- function(sce, celltype.of.interest, microbe.of.interest, sample.of.interest) {
  ci.positive <- sum(counts(sce[, (sce[[celltype.col]] == celltype.of.interest) & (sce[["sample"]] == sample.of.interest)])[microbe.of.interest,] >= min.umis)
  ci.total <- dim(sce[,(sce[[celltype.col]] == celltype.of.interest) & (sce[["sample"]] == sample.of.interest)])[[2]]
  ci.percentage <- ci.positive/ci.total
  return(ci.percentage)
}


Monocyte.Salmonella.exposed <- get.positive.percent(sce, "Monocyte", "Salmonella", "GSM3454529")
nonMonocyte.Salmonella.exposed <- get.positive.percent(sce, "nonMonocyte", "Salmonella", "GSM3454529")
Monocyte.Salmonella.unexposed <- get.positive.percent(sce, "Monocyte", "Salmonella", "GSM3454528")
nonMonocyte.Salmonella.unexposed <- get.positive.percent(sce, "nonMonocyte", "Salmonella", "GSM3454528")


df <- data.frame(celltype=c("Monocyte", "nonMonocyte", "Monocyte", "nonMonocyte"), condition=c("exposed (3' v2)", "exposed (3' v2)", "unexposed (3' v2)", "unexposed (3' v2)"), percentage=c(Monocyte.Salmonella.exposed, nonMonocyte.Salmonella.exposed, Monocyte.Salmonella.unexposed, nonMonocyte.Salmonella.unexposed))
df$celltype <- factor(df$celltype, levels=c("nonMonocyte", "Monocyte"))
df$condition <- factor(df$condition, levels=c("unexposed (3' v2)", "exposed (3' v2)"))
# data <- data.frame(group1=c(celltype.of.interest), group2=c(celltype.comparison), p=c(0.002), y.position=max(nonMonocyte.Salmonella, Monocyte.Salmonella)*1.1)

ggplot(df, aes(fill=celltype, y=percentage, x=condition)) +
    geom_bar(position="dodge", stat="identity") +
    theme_pubr(base_size=6) + # , base_size= 10 to set font size
    # theme(plot.title = element_text(size = 20), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
    ylab("% of cells with >= 1 UMI from Salmonella") +
    xlab(NULL) +
    scale_fill_manual(values=c("#4DBBD5FF", "#E64B35FF")) + # inspired by scale_color_npg
    scale_y_continuous(labels = scales::percent_format(accuracy = .1)) +
    ggtitle("in vitro S. enterica infection (Ben-Moshe2018)")


ggsave(snakemake@output[[1]], width = 2, height = 2, units = "in")
# ggsave("output/plots/figure_1B.pdf", width = 4, height = 6, units = "in")
