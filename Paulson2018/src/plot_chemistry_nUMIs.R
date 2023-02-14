library(scater)
library(scran)
library(ggplot2)
library(ggpubr)
library(ggsci)


counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
# counts <- read.table("output/All/species_PathSeq_All_reads.tsv", sep="\t", header=TRUE, row.names=1)
# output/61_species_PathSeq_metadata.tsv
pdata <- read.table(snakemake@input[[2]], sep="\t", header = TRUE, row.names=1)
# pdata <- read.table("output/All/PathSeq_metadata.tsv", sep="\t", header = TRUE, row.names=1)
row.names(pdata) <- gsub("-", ".", row.names(pdata))
tax.map <- read.table(snakemake@input[[3]], sep="\t", header=TRUE)
# tax.map <- read.table("output/All/tax_id_map_All_PathSeq.tsv", sep="\t", header=TRUE)

column.names <- intersect(unique(colnames(counts)), unique(row.names(pdata)))
# remove samples with no spike-ins

counts <- counts[, column.names]
pdata <- pdata[column.names,]

sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData=pdata)

rownames(sce) <- lapply(rownames(sce), function(x) tax.map[tax.map$tax_id == x, "name"])
microbe.of.interest <- "Human_polyomavirus_5"
celltype.col <- "Tumor"

# get.nreads <- function(sce, sample.of.interest) {
#   nreads <- sum(counts(sce[,(sce[["sample"]] == sample.of.interest)])[microbe.of.interest,])
#   return(nreads)
# }


get.nreads <- function(sce, celltype.of.interest, microbe.of.interest, sample.of.interest) {
  nreads <- sum(counts(sce[, (sce[[celltype.col]] == celltype.of.interest) & (sce[["sample"]] == sample.of.interest)])[microbe.of.interest,])
  return(nreads)
}


Tumor.3p.post <- get.nreads(sce, "Tumor", "Human_polyomavirus_5", "Day615Tumor")
Tumor.3p.pre <- get.nreads(sce, "Tumor", "Human_polyomavirus_5", "PreRxTumor")
Tumor.5p.post <- get.nreads(sce, "Tumor", "Human_polyomavirus_5", "RelapseD565Tumor")
condition <- c("2586 pre", "2586 post", "9245 post")
chemistry <- c("3' v2", "3' v2", "5'")
nreads <- c(Tumor.3p.pre, Tumor.3p.post, Tumor.5p.post)

df <- data.frame(condition, chemistry, nreads)
df$condition <- factor(df$condition, levels=c("2586 pre", "2586 post", "9245 post"))

ggplot(df, aes(fill=chemistry, y=nreads, x=condition)) +
    geom_bar(position="dodge", stat="identity") +
    theme_pubr(base_size=6) + # , base_size= 10 to set font size
    #theme(plot.title = element_text(size = 20), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
    xlab(NULL) +
    ylab("Total number of UMIs to Merkel polyomavirus") +
    scale_fill_manual(values=c("#4DBBD5FF", "#E64B35FF")) + # inspired by scale_color_npg
    ggtitle("Merkel cell carcinoma (Paulson2018)")

# ggsave("output/plots/Fig4C.pdf", width = 5, height = 6, units = "in")
ggsave(snakemake@output[[1]], width = 2, height = 2, units = "in")
