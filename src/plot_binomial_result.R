library(scater)
library(scran)
library(ggplot2)
library(ggpubr)
library(ggforce)

# microbe.counts <- read.table("output/Pt0/GSM3454529_genus_PathSeq_Bacteria_reads.tsv", sep="\t", header=TRUE, row.names=1)
microbe.counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)

# pdata <- read.table("output/Pt0/GSM3454529_genus_PathSeq_Bacteria_metadata.tsv", sep="\t", header=TRUE)
pdata <- read.table(snakemake@input[[2]], sep="\t", header=TRUE)

tax.map <- read.table(snakemake@input[[3]], sep="\t", header=TRUE)

sce <- SingleCellExperiment(assays = list(counts = as.matrix(microbe.counts)), colData=pdata)

rownames(sce) <- lapply(rownames(sce), function(x) tax.map[tax.map$tax_id == x, "name"])

# celltype.col <- "Monocyte"
celltype.col <- snakemake@wildcards[["celltype"]]
# celltype.of.interest <- "Monocyte"
celltype.of.interest <- snakemake@wildcards[["celltype_of_interest"]]
# celltype.comparison <- "nonMonocyte"
celltype.comparison <- snakemake@wildcards[["celltype_comparison"]]
# microbe.of.interest <- "Salmonella"
microbe.of.interest <- snakemake@wildcards[["microbe_of_interest"]]

# now subset the cells by celltypes of interest if applicable
if (celltype.comparison != "all"){
  sce <- sce[, sce[[celltype.col]] %in% c(celltype.of.interest, celltype.comparison)]
}

celltype <- sce[[celltype.col]]

groups <- celltype
# lfc <- 0
lfc <- as.numeric(snakemake@wildcards[["lfc"]])
# pval.type <- "any"
pval.type <- snakemake@wildcards[["pvaltype"]]
# block <- sce[["sample"]]
block <- sce[[snakemake@wildcards[["block"]]]]

binom <- findMarkers(sce, test="binom", groups=groups, pval.type=pval.type, block=block, assay.type="counts")

ci.binom <- binom[[celltype.of.interest]]

# let's create a bar plot of the % of cells with one or more Salmonella UMI
ci.positive <- sum(counts(sce[,groups == celltype.of.interest])[microbe.of.interest,] > 0)
ci.total <- dim(sce[,groups == celltype.of.interest])[[2]]
comparison.positive <- sum(counts(sce[,groups == celltype.comparison])[microbe.of.interest,] > 0)
comparison.total <- dim(sce[,groups == celltype.comparison])[[2]]
ci.percentage <- ci.positive/ci.total
comparison.percentage <- comparison.positive/comparison.total
print(ci.percentage)
print(comparison.percentage)
df <- data.frame(celltype=c(celltype.of.interest, celltype.comparison), percentage=c(ci.percentage, comparison.percentage))
data <- data.frame(group1=c(celltype.of.interest), group2=c(celltype.comparison), p=round(ci.binom[microbe.of.interest, "p.value"], digits=4), y.position=max(ci.percentage, comparison.percentage)*1.1)#

ggbarplot(df, x = "celltype", y = "percentage",
   fill = "celltype", color = "celltype") +
   scale_fill_manual(values=c("#e74c3c", "#3498db")) +
   ylab(paste0("Percentage of cells with UMIs mapped to ", microbe.of.interest)) +
   xlab("Celltype") +
   ggtitle(paste0("Patient ", snakemake@wildcards[["patient"]])) +
   theme_pubr(legend="none") +
   stat_pvalue_manual(data, "binomial test p-value = {p}") +
   scale_y_continuous(labels = scales::percent_format(accuracy = .1)) +
   font("xy.text", size=14) +
   font("xylab", size=14) +
   theme(plot.title = element_text(hjust = 0.5)) +
   font("title", size=18)

ggsave(snakemake@output[[1]])
