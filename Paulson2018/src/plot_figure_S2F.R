library(scater)
library(scran)
library(ggplot2)
library(ggpubr)
library(ggsci)

# counts <- read.table("output/9245-3/RelapseD565Tumor/genus_PathSeq_All_reads.tsv", sep="\t", header=TRUE, row.names=1)
counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
# pdata <- read.table("output/9245-3/PathSeq_metadata.tsv", sep="\t", header = TRUE, row.names=1)
pdata <- read.table(snakemake@input[[2]], sep="\t", header = TRUE, row.names=1)

row.names(pdata) <- gsub("-", ".", row.names(pdata))
# tax.map <- read.table("output/9245-3/tax_id_map_All_PathSeq.tsv", sep="\t", header=TRUE)
tax.map <- read.table(snakemake@input[[3]], sep="\t", header=TRUE)

column.names <- intersect(unique(colnames(counts)), unique(row.names(pdata)))
# remove samples with no spike-ins

counts <- counts[, column.names]
pdata <- pdata[column.names,]

sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData=pdata)

rownames(sce) <- lapply(rownames(sce), function(x) tax.map[tax.map$tax_id == x, "name"])
print(sce)
microbe.of.interest <- "Human_polyomavirus_5"
min.umis <- 2 #strtoi(snakemake@wildcards[["min_umis"]])

get.positive.percent <- function(sce, celltype.of.interest) {
  ci.positive <- sum(counts(sce[,(sce[["Tumor"]] == celltype.of.interest)])[microbe.of.interest,] >= min.umis)
  ci.total <- dim(sce[,(sce[["Tumor"]] == celltype.of.interest)])[[2]]
  ci.percentage <- ci.positive/ci.total
  return(ci.percentage)
}



Tumor.5p.post <- get.positive.percent(sce, "Tumor")
non.Tumor.5p.post <- get.positive.percent(sce, "nonTumor")
condition <- c("9245 post (5')", "9245 post (5')")
celltype <- c("Tumor", "nonTumor")
positive.percent <- c(Tumor.5p.post, non.Tumor.5p.post)

df <- data.frame(condition, celltype, positive.percent)
print(df)
df$condition <- factor(df$condition, levels=c("9245 post (5')"))
df$celltype <- factor(df$celltype, levels=c("nonTumor", "Tumor"))

# print(df)
ggplot(df, aes(fill=celltype, y=positive.percent, x=celltype)) +
    geom_bar(position="dodge", stat="identity") +
    theme_pubr(base_size=6, legend="none") + # , base_size= 10 to set font size
    # theme(plot.title = element_text(size = 20), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
    xlab(NULL) +
    ylab("% cells with >= 2 UMIs from Merkel polyomavirus") +
    scale_fill_manual(values=c("#4DBBD5FF", "#E64B35FF")) + # inspired by scale_color_npg
    scale_y_continuous(labels = scales::percent_format(accuracy = .1), limits=c(0, .50)) +
    ggtitle("9245 post (5')")

# ggsave("output/plots/Fig4C.pdf", width = 5, height = 6, units = "in")
ggsave(snakemake@output[[1]], width = 2, height = 2, units = "in")
