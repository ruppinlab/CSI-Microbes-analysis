library(scater)
library(scran)
library(ggplot2)
library(ggpubr)
library(ggsci)


counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
# counts <- read.table("output/Pt0/GSM3454529/genus_PathSeq_Bacteria_reads.tsv", sep="\t", header=TRUE, row.names=1)
# output/61_species_PathSeq_metadata.tsv
pdata <- read.table(snakemake@input[[2]], sep="\t", header = TRUE, row.names=1)
row.names(pdata) <- gsub("-", ".", row.names(pdata))
# pdata <- read.table("output/Pt0/GSM3454529/genus_PathSeq_Bacteria_metadata.tsv", sep="\t", header = TRUE, row.names=1)
tax.map <- read.table(snakemake@input[[3]], sep="\t", header=TRUE)
# tax.map <- read.table("output/Pt0/tax_id_map_Bacteria_PathSeq.tsv", sep="\t", header=TRUE)
sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData=pdata)

rownames(sce) <- lapply(rownames(sce), function(x) tax.map[tax.map$tax_id == x, "name"])

celltype.col <- "Tumor"
# celltype.col <- snakemake@wildcards[["celltype"]]
celltype.of.interest <- "Tumor"
# celltype.of.interest <- snakemake@wildcards[["celltype_of_interest"]]
celltype.comparison <- "nonTumor"
# celltype.comparison <- snakemake@wildcards[["celltype_comparison"]]
microbe.of.interest <- "Fusobacterium_nucleatum"
# microbe.of.interest <- snakemake@wildcards[["microbe_of_interest"]]

# now subset the cells by celltypes of interest if applicable
if (celltype.comparison != "all"){
  sce <- sce[, sce[[celltype.col]] %in% c(celltype.of.interest, celltype.comparison)]
}

celltype <- sce[[celltype.col]]
#
# counts <- t(counts(sce))
# colnames(logcounts) <- lapply(colnames(logcounts), function(x) tax.map[tax.map$tax_id == x, "name"])
# df <- cbind(as.data.frame(t(counts(sce))), as.data.frame(colData(sce)))

ci.positive <- sum(counts(sce[,celltype == celltype.of.interest])[microbe.of.interest,] > 0)
ci.total <- dim(sce[,celltype == celltype.of.interest])[[2]]
comparison.positive <- sum(counts(sce[,celltype == celltype.comparison])[microbe.of.interest,] > 0)
comparison.total <- dim(sce[,celltype == celltype.comparison])[[2]]
ci.percentage <- ci.positive/ci.total
comparison.percentage <- comparison.positive/comparison.total
#print(ci.percentage)
#print(comparison.percentage)
df <- data.frame(celltype=c(celltype.comparison, celltype.of.interest), percentage=c(comparison.percentage, ci.percentage))
df$celltype <- factor(df$celltype, levels=c("nonTumor", "Tumor"))
#data <- data.frame(group1=c(celltype.of.interest), group2=c(celltype.comparison), p=round(ci.binom[microbe.of.interest, "p.value"], digits=4), y.position=max(ci.percentage, comparison.percentage)*1.1)#



fig.1d <- ggplot(df, aes(celltype, percentage, fill=celltype)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  theme_pubr(base_size=8) + # , base_size= 10 to set font size
  ylab("% of cells with >= 1 UMI") +
  xlab(NULL) +
  # scale_color_npg() +
  scale_fill_manual(values=c("#4DBBD599", "#E64B3599")) + # inspired by scale_color_npg
  scale_y_continuous(labels = scales::percent_format(accuracy = .1)) +
  ggtitle("Fusobacterium nucleatum")

#
# ggbarplot(df, x = "celltype", y = "percentage",
#    fill = "celltype", color = "celltype") +
#    scale_fill_manual(values=c("#e74c3c", "#3498db")) +
#    ylab("Percentage of cells") +
#    xlab("Celltype") +
#    ggtitle("Cells with UMIs mapped to Salmonella") +
#    theme_pubr(legend="none") +
#    stat_pvalue_manual(data, "binomial test p-value = {p}") +
#    scale_y_continuous(labels = scales::percent_format(accuracy = 1))

ggsave(snakemake@output[[1]], width = 2, height = 3, units = "in")
