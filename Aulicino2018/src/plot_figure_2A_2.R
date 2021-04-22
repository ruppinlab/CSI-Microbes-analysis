library(ggplot2)
library(ggpubr)
library(ggsci)
library(scater)
library(scran)


source(snakemake@params[["spike_functions"]])

# sce <- generate_sce(microbe.file="output/Pt0/class_PathSeq_Bacteria_reads.tsv",
# spikein.file="output/star/Pt0/spike_ins.tsv", pdata.file="output/Pt0/class_PathSeq_Bacteria_metadata.tsv", celltype.col="infected")
sce <- generate_sce(microbe.file=snakemake@input[[1]], spikein.file=snakemake@input[[2]],
  pdata.file=snakemake@input[[3]], celltype.col="infected")

# sce <- sce[, sce$infected %in% c("infected", "bystander")]


tax.map <- read.table(snakemake@input[[4]], sep="\t", header=TRUE)
# tax.map <- read.table("output/Pt0/tax_id_map_Bacteria_PathSeq.tsv", sep="\t", header=TRUE)
#species.tax.map <- tax.map[tax.map$taxa_level == snakemake@wildcards[["tax_level"]],]

logcounts <- t(logcounts(sce))
colnames(logcounts) <- lapply(colnames(logcounts), function(x) tax.map[tax.map$tax_id == x, "name"])
df <- cbind(as.data.frame(logcounts), as.data.frame(colData(sce)))

#print(df)

# df <- cbind(pdata, counts)
# df$logcounts <- log2(df$counts + 1)
df$infected <- factor(df$infected, levels=c("uninfected", "bystander", "infected"))
levels(df$infected) <- c("Control", "Bystander", "Infected")

# print(pal_npg("nrc")(3))
# print(pal_npg("nrc", alpha=0.6)(3))

# rename rownames from tax_id to names
fig.1a <- ggplot(df, aes(x=infected, y=Salmonella_enterica)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(aes(color=infected), alpha=0.6, size=2) +
  theme_pubr(legend="none", base_size=9) + # , base_size= 10 to set font size
  ylab("log2(spike-in normalized abundance + 1)") +
  xlab(NULL) +
  # scale_color_npg() +
  scale_color_manual(values=c("#00A087FF", "#4DBBD5FF","#E64B35FF")) + # inspired by scale_color_npg
  ggtitle("Salmonella enterica abundance")

print(ggplot_build(fig.1a)$data)

# ggtitle("Reads mapped to \nSalmonella genome") +

ggsave(snakemake@output[[1]], width = 2.8, height = 3, units = "in")
