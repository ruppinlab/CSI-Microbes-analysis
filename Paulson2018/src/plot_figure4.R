#library(EnhancedVolcano)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggrepel)

# taken from https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html

my_data <- lapply(snakemake@input, read.table)
#my_data <- lapply(input, read.table)
print(dim(my_data))

df <- do.call(rbind, my_data)


# EnhancedVolcano(df, lab=df$taxa, x = "summary.logFC", y = "p.value",
#   FCcutoff = .5, pCutoff = .05, selectLab = microbes.of.interest, subtitle=NULL,
#   caption=NULL, title="Infected versus Bystander", xlim=c(-2,3),
#   ylim = c(0, -log10(10e-13))
# )

df$diffexpressed <- "NO"
df$diffexpressed[df$summary.logFC > 0.5 & df$p.value < 0.05] <- "UP"
df$diffexpressed[df$summary.logFC < -0.5 & df$p.value < 0.05] <- "DOWN"

df$delabel <- NA
df$delabel[df$diffexpressed != "NO"] <- df$taxa[df$diffexpressed != "NO"]

mycolors <- c("red", "red", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")

ggplot(df, aes(x=summary.logFC, y=-log10(p.value), color=diffexpressed, label=delabel)) +
  geom_point(alpha=0.6, size=2) +
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype="longdash") + #
  geom_hline(yintercept=-log10(0.05), col="black", linetype="longdash") +
  geom_text_repel(color="black", size=3) +
  scale_colour_manual(values = mycolors) +
  theme_pubr(legend="none", base_size=8) +
  ggtitle("Tumor versus non-tumor cells") +
  ylab("-log10(P)") +
  xlab("log2 fold change")

  # ylim(0, 12.5) +
  # xlim(-2, 3) +

ggsave(snakemake@output[[1]], width=3, height=3, units = "in")
