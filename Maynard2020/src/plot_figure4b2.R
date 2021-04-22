library(ggplot2)
library(ggpubr)
library(ggsci)
library(scater)
library(scran)
library(stringr)

# source("../src/spike_normalization_functions.R")
source(snakemake@params[["spike_functions"]])

celltype.col <- "celltype1"
# celltype.col <- snakemake@wildcards[["celltype"]]
celltype.of.interest <- "Tumor"
# celltype.of.interest <- snakemake@wildcards[["celltype_of_interest"]]
celltype.comparison <- "immune"
# celltype.comparison <- snakemake@wildcards[["celltype_comparison"]]

microbe.of.interest <- "Cutibacterium_acnes"
# microbe.of.interest <- snakemake@wildcards[["microbe_of_interest"]]

# sce <- generate_sce(microbe.file="output/TH238/species_PathSeq_Bacteria_reads.tsv",
# spikein.file="output/star/TH238/spike_ins.tsv", pdata.file="output/TH238/species_PathSeq_Bacteria_metadata.tsv", celltype.col="celltype1")
sce <- generate_sce(microbe.file=snakemake@input[[1]], spikein.file=snakemake@input[[2]],
  pdata.file=snakemake@input[[3]], celltype.col=celltype.col)

# sce <- sce[, sce$infected %in% c("infected", "bystander")]


tax.map <- read.table(snakemake@input[[4]], sep="\t", header=TRUE)
# tax.map <- read.table("output/TH238/tax_id_map_Bacteria_PathSeq.tsv", sep="\t", header=TRUE)
#species.tax.map <- tax.map[tax.map$taxa_level == snakemake@wildcards[["tax_level"]],]

sce <- sce[, sce[[celltype.col]] %in% c(celltype.of.interest, celltype.comparison)]

logcounts <- t(logcounts(sce))
colnames(logcounts) <- lapply(colnames(logcounts), function(x) tax.map[tax.map$tax_id == x, "name"])
df <- cbind(as.data.frame(logcounts), as.data.frame(colData(sce)))

#print(df)


# df <- cbind(pdata, counts)
# df$logcounts <- log2(df$counts + 1)
df[[celltype.col]] <- factor(df[[celltype.col]], levels=c(celltype.comparison, celltype.of.interest))
levels(df[[celltype.col]]) <- c(str_to_title(celltype.comparison), str_to_title(celltype.of.interest))

# rename rownames from tax_id to names
fig.1a <- ggplot(df, aes_string(x=celltype.col, y=microbe.of.interest)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(aes_string(color=celltype.col), alpha=0.6, size=2) +
  theme_pubr(legend="none", base_size=8) + # , base_size= 10 to set font size
  ylab("log2(spike-in normalized abundance + 1)") +
  xlab(NULL) +
  # scale_color_npg() +
  scale_color_manual(values=c("#4DBBD5FF","#E64B35FF")) + # inspired by scale_color_npg
  ggtitle("Cutibacterium acnes") +

ggsave(snakemake@output[[1]], width = 2, height = 3, units = "in")
