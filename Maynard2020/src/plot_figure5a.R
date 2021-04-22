library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggrepel)

# taken from https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html

# df <- read.table("output/results/TH236/hFDR-wilcox-celltype1-Tumor-immune-PathSeq-spike-Bacteria-any-.5-plate-up-0.5.tsv")
df <- read.table(snakemake@input[[1]], header=TRUE)

# EnhancedVolcano(df, lab=df$taxa, x = "summary.logFC", y = "p.value",
#   FCcutoff = .5, pCutoff = .05, selectLab = microbes.of.interest, subtitle=NULL,
#   caption=NULL, title="Infected versus Bystander", xlim=c(-2,3),
#   ylim = c(0, -log10(10e-13))
# )

df$diffexpressed <- "NO"
df$diffexpressed[df$summary.AUC > 0.55 & df$p.value < 10e-6] <- "UP"
df$diffexpressed[df$summary.AUC < 0.45 & df$p.value < 10e-6] <- "DOWN"

df$delabel <- NA
df$delabel[df$diffexpressed != "NO"] <- rownames(df)[df$diffexpressed != "NO"]

mycolors <- c("red", "red", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")

ggplot(df, aes(x=summary.AUC, y=-log10(p.value), color=diffexpressed, label=delabel)) +
  geom_point(alpha=0.6, size=2, show.legend = FALSE) +
  geom_vline(xintercept=c(0.45, 0.55), col="black", linetype="longdash") +
  geom_hline(yintercept=-log10(10e-6), col="black", linetype="longdash") +
  geom_text_repel(color="black", size=3) +
  scale_colour_manual(values = mycolors) +
  theme_light() +
  ggtitle("infected versus uninfected tumor cells") +
  ylab("-log10(P)") +
  xlab("AUC") +
  ylim(0, 75)

#   theme_pubr(legend="none", base_size=8) +
  # ylim(0, 18) +
  # xlim(-2, 3) +

ggsave(snakemake@output[[1]], width=4, height=4, units = "in")
