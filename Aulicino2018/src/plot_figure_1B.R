library(ggplot2)
library(ggpubr)
library(ggsci)
library(scater)
library(scran)


source(snakemake@params[["spike_functions"]])
# source("../src/spike_normalization_functions.R")

# sce <- generate_sce(microbe.file="output/Pt0/genus_PathSeq_All_reads.tsv",
# spikein.file="output/star/Pt0/spike_ins.tsv", pdata.file="output/Pt0/genus_PathSeq_All_metadata.tsv", celltype.col="infected")
sce <- generate_sce(microbe.file=snakemake@input[[1]], spikein.file=snakemake@input[[2]],
  pdata.file=snakemake@input[[3]], celltype.col="infected")

# sce <- sce[, sce$infected %in% c("infected", "bystander")]


tax.map <- read.table(snakemake@input[[4]], sep="\t", header=TRUE)
# tax.map <- read.table("output/Pt0/tax_id_map_Bacteria_PathSeq.tsv", sep="\t", header=TRUE)
#species.tax.map <- tax.map[tax.map$taxa_level == snakemake@wildcards[["tax_level"]],]

logcounts <- t(logcounts(sce))
colnames(logcounts) <- lapply(colnames(logcounts), function(x) tax.map[tax.map$tax_id == x, "name"])
# create data.frame
# df <- cbind(as.data.frame(logcounts), as.data.frame(colData(sce)))
df <- data.frame(infected=sce$infected, Salmonella=logcounts[, "Salmonella"])
df$infected <- factor(df$infected, levels=c("uninfected", "bystander", "infected"))
levels(df$infected) <- c("Control", "Bystander", "Infected")

max.Salmonella.infected.value <- max(df[df$infected == "Infected", "Salmonella"])
max.Salmonella.bystander.value <- max(df[df$infected == "Bystander", "Salmonella"])

data <- data.frame(group1=c("Infected", "Infected", "Bystander"), group2=c("Bystander", "Control", "Control"), p=c("***", "***", "***"), y.position=c(max.Salmonella.infected.value*1.1, max.Salmonella.infected.value*1.05, max.Salmonella.bystander.value*1.05))

# print(pal_npg("nrc")(3))
# print(pal_npg("nrc", alpha=0.6)(3))

# rename rownames from tax_id to names
# fig.1a <- ggplot(df, aes(x=infected, y=Salmonella)) +
#   geom_boxplot(outlier.alpha = 0) +
#   geom_jitter(aes(color=infected), alpha=0.6, size=2) +
#   theme_pubr(legend="none", base_size=9) + # , base_size= 10 to set font size
#   ylab("log2(spike-in normalized abundance + 1)") +
#   xlab(NULL) +
#   # scale_color_npg() +
#   scale_color_manual(values=c("#00A087FF", "#4DBBD5FF","#E64B35FF")) + # inspired by scale_color_npg
#   ggtitle("Salmonella abundance")

fig.1a <- ggboxplot(df, x = "infected", y = "Salmonella", color = "infected", add = "jitter") +
  # geom_jitter(aes(color=infected), alpha=0.6, size=2) +
  theme_pubr(legend="none", base_size=6) + # , base_size= 10 to set font size
  ylab("log2(normalized Salmonella reads + 1)") +
  xlab(NULL) +
  # theme(plot.title = element_text(size = 20), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  stat_pvalue_manual(data, label="{p}", tip.length = 0, size=2) +
  scale_color_manual(values=c("#00A087FF", "#4DBBD5FF","#E64B35FF")) + # inspired by scale_color_npg
  ggtitle("in vitro S. enterica infection (Aulicino2018)")

# print(ggplot_build(fig.1a)$data)

# ggtitle("Reads mapped to \nSalmonella genome") +

ggsave(snakemake@output[[1]], width = 2, height = 2, units = "in")
