library(scater)
library(scran)
library(ggplot2)
library(ggpubr)
library(ggsci)


counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
# counts <- read.table("output/C166/C166_T_1_1_0_c1_v3/genus_PathSeq_All_reads.tsv", sep="\t", header=TRUE, row.names=1)
# output/61_species_PathSeq_metadata.tsv
pdata <- read.table(snakemake@input[[2]], sep="\t", header = TRUE, row.names=1)
# pdata <- read.table("output/C163/PathSeq_metadata.tsv", sep="\t", header = TRUE, row.names=1)
row.names(pdata) <- gsub("-", ".", row.names(pdata))

tax.map <- read.table(snakemake@input[[3]], sep="\t", header=TRUE)
# tax.map <- read.table("output/C163/tax_id_map_All_PathSeq.tsv", sep="\t", header=TRUE)

column.names <- intersect(unique(colnames(counts)), unique(row.names(pdata)))
# remove samples with no spike-ins

counts <- counts[, column.names]
pdata <- pdata[column.names,]

sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData=pdata)

rownames(sce) <- lapply(rownames(sce), function(x) tax.map[tax.map$tax_id == x, "name"])

# microbe.of.interest <- "Fusobacterium"
min.umis <- 2 #strtoi(snakemake@wildcards[["min_umis"]])

get.positive.percent <- function(sce, celltype.of.interest, microbe.of.interest) {
  ci.positive <- sum(counts(sce[,(sce[["celltype1"]] == celltype.of.interest)])[microbe.of.interest,] >= min.umis)
  ci.total <- dim(sce[,(sce[["celltype1"]] == celltype.of.interest)])[[2]]
  ci.percentage <- ci.positive/ci.total
  return(ci.percentage)
}


Myeloid.Fusobacterium <- get.positive.percent(sce, "Myeloid", "Fusobacterium")
Tcell.Fusobacterium <- get.positive.percent(sce, "TNKILC", "Fusobacterium")
Bcell.Fusobacterium <- get.positive.percent(sce, "B", "Fusobacterium")
Epithelial.Fusobacterium <- get.positive.percent(sce, "Epithelial", "Fusobacterium")
Stromal.Fusobacterium <- get.positive.percent(sce, "Stromal", "Fusobacterium")


condition <- c(rep("Fusobacterium", 5))
celltype <- c(rep(c("Myeloid", "TNKILC", "B", "Tumor", "Stromal"), 3))
positive.percent <- c(Myeloid.Fusobacterium, Tcell.Fusobacterium, Bcell.Fusobacterium, Epithelial.Fusobacterium, Stromal.Fusobacterium)

df <- data.frame(condition, celltype, positive.percent)

df$condition <- factor(df$condition, levels=c("Fusobacterium"))
df$celltype <- factor(df$celltype, levels=c("Myeloid", "Stromal", "Tumor", "TNKILC", "B"))

ggplot(df, aes(fill=celltype, y=positive.percent, x=condition)) +
    geom_bar(position="dodge", stat="identity") +
    theme_pubr(base_size=6) + # , base_size= 10 to set font size
    ylab("% of cells infected") +
    theme(legend.key.size = unit(.1, 'in')) +
    scale_color_npg() +
    scale_y_continuous(labels = scales::percent_format(accuracy = .1)) +
    ggtitle("C169")

# scale_fill_manual() + # inspired by scale_color_npg - values=c("#00A087FF", "#4DBBD5FF","#E64B35FF")

ggsave(snakemake@output[[1]], width = 2.5, height = 2, units = "in")
