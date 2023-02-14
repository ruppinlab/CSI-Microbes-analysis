library(ggplot2)
library(ggpubr)
library(ggsci)
library(scater)
library(scran)

# source("../src/spike_normalization_functions.R")
source(snakemake@params[["spike_functions"]])

sce <- generate_sce(microbe.file="output/P1/genus_PathSeq_All_reads.tsv",
spikein.file="output/star/P1/spike_ins.tsv", pdata.file="output/P1/genus_PathSeq_All_metadata.tsv", celltype.col="Cell_Type_and_Condition")
sce <- generate_sce(microbe.file=snakemake@input[[1]], spikein.file=snakemake@input[[2]],
  pdata.file=snakemake@input[[3]], celltype.col="Cell_Type_and_Condition")

print(sce)

sce <- sce[, sce$Cell_Type_and_Condition %in% c("Jurkat-uninfected", "HCT116-uninfected", "Jurkat-heat-killed", "HCT116-heat-killed", "HCT116-infected", "Jurkat-infected", "ERCC only")]

print(sce)

tax.map <- read.table(snakemake@input[[4]], sep="\t", header=TRUE)
# tax.map <- read.table("output/P1/tax_id_map_All_PathSeq.tsv", sep="\t", header=TRUE)


logcounts <- t(logcounts(sce))
colnames(logcounts) <- lapply(colnames(logcounts), function(x) tax.map[tax.map$tax_id == x, "name"])

df <- data.frame(Cell_Type_and_Condition=sce$Cell_Type_and_Condition, Condition=sce$Condition, Cell_Type=sce$Cell_Type, Fusobacterium=logcounts[, "Fusobacterium"])
print(df)
df$Condition <- replace(df$Condition, df$Condition == "ERCC only", "Empty Wells")
df$Condition <- replace(df$Condition, df$Condition == "Uninfected", "Unexposed")
df$Condition <- replace(df$Condition, df$Condition == "Heat Killed", "Heat-Killed Fn")
df$Condition <- replace(df$Condition, df$Condition == "Infected", "Live Fn")

df$Cell_Type_and_Condition <- replace(df$Cell_Type_and_Condition, df$Cell_Type_and_Condition == "ERCC only", "Empty Wells")
df$Cell_Type <- replace(df$Cell_Type, df$Cell_Type == "ERCC only", "Empty Wells")

print(df)

df$Cell_Type_and_Condition <- factor(df$Cell_Type_and_Condition, levels=c("Empty Wells", "Jurkat-uninfected", "HCT116-uninfected", "Jurkat-heat-killed", "HCT116-heat-killed", "Jurkat-infected", "HCT116-infected"))

df$Cell_Type <- factor(df$Cell_Type, levels=c("Empty Wells", "Jurkat", "HCT116"))
df$Condition <- factor(df$Condition, levels=c("Empty Wells", "Unexposed", "Heat-Killed Fn", "Live Fn"))
# three cell-types - Empty Wells, Jurkat T and HCT116
# three conditions - uninfected, heat-killed and live exposed

#max.Fusobacterium.infected.value <- max(df[df$Condition == "Infected", "Fusobacterium"])
#max.Fusobacterium.HK.value <- max(df[df$Condition == "Heat-Killed", "Fusobacterium"])


#data <- data.frame(group1=c("Infected", "Infected", "Infected"), group2=c("Heat-Killed", "Uninfected", "Empty Wells"), p=c(.00676628201353454, 0.00285014808429544, 0.000233521033034032), y.position=c(max.Fusobacterium.infected.value*1.05, max.Fusobacterium.infected.value*1.1, max.Fusobacterium.infected.value*1.15))


fig.1a <- ggboxplot(df, x = "Condition", y = "Fusobacterium", color = "Cell_Type", add = "jitter") +
  # geom_jitter(aes(color=infected), alpha=0.6, size=2) +
  theme_pubr(base_size=6) + # , base_size= 10 to set font size
  ylab("log2(normalized Fusobacterium reads + 1)") +
  # theme(plot.title = element_text(size = 20), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  xlab(NULL) +
  # stat_pvalue_manual(data, label="{p}", tip.length = 0, size=2) +
  scale_color_manual(values=c("#00a087", "#4DBBD5FF","#E64B35FF")) + # inspired by scale_color_npg
  ggtitle("plate-based scRNA-seq (Robinson2023)")


# ggsave("output/plots/Figure3A.pdf", width = 7, height = 6, units = "in")
ggsave(snakemake@output[[1]], width = 3, height = 2, units = "in")
