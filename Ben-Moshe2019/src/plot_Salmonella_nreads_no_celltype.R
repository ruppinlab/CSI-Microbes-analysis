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

# celltype <- sce[[celltype.col]]
#
# counts <- t(counts(sce))
# colnames(logcounts) <- lapply(colnames(logcounts), function(x) tax.map[tax.map$tax_id == x, "name"])
# df <- cbind(as.data.frame(t(counts(sce))), as.data.frame(colData(sce)))


get.nreads <- function(sce, microbe.of.interest, sample.of.interest) {
  n.reads <- sum(counts(sce[, (sce[["sample"]] == sample.of.interest)])[microbe.of.interest,])
  return(log2(n.reads+1))
}

Salmonella.exposed <- get.nreads(sce, "Salmonella", "GSM3454529")
Salmonella.unexposed <- get.nreads(sce, "Salmonella", "GSM3454528")
print(Salmonella.exposed)
print(Salmonella.unexposed)
df <- data.frame(condition=c("unexposed", "exposed"), nreads=c(Salmonella.unexposed, Salmonella.exposed))
df$condition <- factor(df$condition, levels=c("unexposed", "exposed"))
# data <- data.frame(group1=c(celltype.of.interest), group2=c(celltype.comparison), p=c(0.002), y.position=max(nonMonocyte.Salmonella, Monocyte.Salmonella)*1.1)

ggplot(df, aes(fill=condition, y=nreads, x=condition)) +
    geom_bar(position="dodge", stat="identity", show.legend = FALSE) +
    theme_pubr(base_size=6) + # , base_size= 10 to set font size
    # theme(plot.title = element_text(size = 20), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
    ylab("number of UMIs from Salmonella") +
    xlab(NULL) +
    scale_fill_manual(values=c("#4DBBD5FF", "#E64B35FF")) + # inspired by scale_color_npg
    ggtitle("in vitro S. enterica infection (Ben-Moshe2018)")


ggsave(snakemake@output[[1]], width = 2, height = 2, units = "in")
# ggsave("output/plots/figure_1B.pdf", width = 4, height = 6, units = "in")
