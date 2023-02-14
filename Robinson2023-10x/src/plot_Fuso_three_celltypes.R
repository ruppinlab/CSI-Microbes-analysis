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

column.names <- intersect(unique(colnames(counts)), unique(row.names(pdata)))
# remove samples with no spike-ins

counts <- counts[, column.names]
pdata <- pdata[column.names,]

sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData=pdata)

rownames(sce) <- lapply(rownames(sce), function(x) tax.map[tax.map$tax_id == x, "name"])

celltype.col <- "celltype1"
# celltype.col <- snakemake@wildcards[["celltype"]]
celltype1 <- snakemake@wildcards[["celltype1"]]
celltype2 <- snakemake@wildcards[["celltype2"]]
celltype3 <- snakemake@wildcards[["celltype3"]]
microbe.of.interest <- "Fusobacterium"
min.umis <- strtoi(snakemake@wildcards[["min_umis"]])

# now subset the cells by celltypes of interest if applicable

sce <- sce[, sce[[celltype.col]] %in% c(celltype1, celltype2, celltype3)]


celltype <- sce[[celltype.col]]

c1.positive <- sum(counts(sce[,celltype == celltype1])[microbe.of.interest,] >= min.umis)
c1.total <- dim(sce[,celltype == celltype1])[[2]]
c2.positive <- sum(counts(sce[,celltype == celltype2])[microbe.of.interest,] >= min.umis)
c2.total <- dim(sce[,celltype == celltype2])[[2]]
c3.positive <- sum(counts(sce[,celltype == celltype3])[microbe.of.interest,] >= min.umis)
c3.total <- dim(sce[,celltype == celltype3])[[2]]

c1.percentage <- c1.positive/c1.total
c2.percentage <- c2.positive/c2.total
c3.percentage <- c3.positive/c3.total

#print(ci.percentage)
#print(comparison.percentage)
df <- data.frame(celltype=c(celltype1, celltype2, celltype3), percentage=c(c1.percentage, c2.percentage, c3.percentage))
df$celltype <- factor(df$celltype, levels=c(celltype1, celltype2, celltype3))
#data <- data.frame(group1=c(celltype.of.interest), group2=c(celltype.comparison), p=round(ci.binom[microbe.of.interest, "p.value"], digits=4), y.position=max(ci.percentage, comparison.percentage)*1.1)#

# stat.test <- tibble::tribble(
#   ~group1, ~group2, ~p.adj,
#   "Jurkat", "HCT116", 3.1115424659124e-11,
# )

fig.1d <- ggplot(df, aes(celltype, percentage, fill=celltype)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  theme_pubr(base_size=8) + # , base_size= 10 to set font size
  ylab(paste0("% of cells >= ", min.umis, " UMI")) +
  xlab(NULL) +
  scale_color_npg() +
  # scale_fill_manual(values=c("#4DBBD599", "#E64B3599")) + # inspired by scale_color_npg
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

ggsave(snakemake@output[[1]], width = 3, height = 3, units = "in")
