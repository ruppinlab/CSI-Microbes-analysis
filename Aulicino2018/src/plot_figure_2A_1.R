#library(EnhancedVolcano)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggrepel)

# taken from https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html

# input <- c("output/Pt0/wilcox-infected-infected-bystander-species-PathSeq-spike-All-any-0.5-plate-both-0.tsv", "output/Pt0/wilcox-infected-infected-bystander-genus-PathSeq-spike-All-any-0.5-plate-both-0.tsv", "output/Pt0/wilcox-infected-infected-bystander-family-PathSeq-spike-All-any-0.5-plate-both-0.tsv", "output/Pt0/wilcox-infected-infected-bystander-order-PathSeq-spike-All-any-0.5-plate-both-0.tsv", "output/Pt0/wilcox-infected-infected-bystander-class-PathSeq-spike-All-any-0.5-plate-both-0.tsv")
my_data <- lapply(snakemake@input, read.table, header=TRUE)
#my_data <- lapply(input, read.table)


df <- do.call(rbind, my_data)

df$diffexpressed <- "NO"
df$diffexpressed[df$summary.AUC > 0.5 & df$p.value < 0.05] <- "UP"
df$diffexpressed[df$summary.AUC < 0.5 & df$p.value < 0.05] <- "DOWN"

df$delabel <- NA
df$delabel[df$diffexpressed != "NO"] <- df$taxa[df$diffexpressed != "NO"]

mycolors <- c("red", "red", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")

ggplot(df, aes(x=summary.AUC, y=-log10(p.value), color=diffexpressed, label=delabel)) +
  geom_point(alpha=0.6, size=2, show.legend = FALSE) +
  geom_vline(xintercept=c(0.5), col="black", linetype="longdash") + #
  geom_hline(yintercept=-log10(0.05), col="black", linetype="longdash") +
  geom_text_repel(color="black", size=3) +
  scale_colour_manual(values = mycolors) +
  ylim(0, 11) +
  xlim(.35, .8) +
  theme_light(base_size=8) +
  ggtitle("Infected versus Bystander cells") +
  ylab("-log10(P)") +
  xlab("AUC")

#   xlim(-2, 3) +

ggsave(snakemake@output[[1]], width=3, height=3, units = "in")
