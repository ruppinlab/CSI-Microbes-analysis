#library(EnhancedVolcano)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggrepel)
library(scales)


# input <- c("output/Pt0/GSM3454529/fisher-exact-Monocyte-Monocyte-nonMonocyte-phylum-PathSeq-Bacteria.tsv", "output/Pt0/GSM3454529/fisher-exact-Monocyte-Monocyte-nonMonocyte-class-PathSeq-Bacteria.tsv", "output/Pt0/GSM3454529/fisher-exact-Monocyte-Monocyte-nonMonocyte-order-PathSeq-Bacteria.tsv", "output/Pt0/GSM3454529/fisher-exact-Monocyte-Monocyte-nonMonocyte-family-PathSeq-Bacteria.tsv", "output/Pt0/GSM3454529/fisher-exact-Monocyte-Monocyte-nonMonocyte-genus-PathSeq-Bacteria.tsv", "output/Pt0/GSM3454529/fisher-exact-Monocyte-Monocyte-nonMonocyte-species-PathSeq-Bacteria.tsv")
# print(snakemake@input)
# print(unlist(snakemake@input))
# print(unlist(snakemake@input)[c(2,3,4,5,6)])
# print(snakemake@input[[2:]])

my_data <- lapply(unlist(snakemake@input)[c(2,3,4,5,6)], read.table)
#my_data <- lapply(input, read.table)


df <- do.call(rbind, my_data)
#print(df)

hFDR.df <- read.table(snakemake@input[[1]], header=TRUE)
print(hFDR.df)

# taken from - https://stackoverflow.com/questions/47085514/simple-way-to-visualise-odds-ratios-in-r
# scales::pseudo_log_trans taken from https://stackoverflow.com/questions/40219639/how-to-deal-with-zero-in-log-plot

# only keep taxa kept by hFDR filtering
df <- df[(df$taxa %in% hFDR.df$name),]

# print(df[order(df$p.value),])
df$taxa <- factor(df$taxa, levels=df[order(-df$p.value),]$taxa)
print(df[order(-df$p.value),])
# levels(df$taxa) <- taxa

ggplot(df, aes(x = summary.odds.ratio, y = taxa)) +
  geom_vline(aes(xintercept = 1), size = .25, linetype="dashed") +
  geom_errorbarh(aes(xmax = or.ci.high, xmin=or.ci.low), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange") +
  theme_pubr(base_size=8) +
  ggtitle("SC028 tumor versus non-tumor cells") +
  ylab(NULL) +
  xlab("odds ratio") +
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10))




# df$diffexpressed <- "NO"
# df$diffexpressed[df$summary.logFC > 0.5 & df$p.value < 0.05] <- "UP"
# df$diffexpressed[df$summary.logFC < -0.5 & df$p.value < 0.05] <- "DOWN"
#
# df$delabel <- NA
# df$delabel[df$diffexpressed != "NO"] <- df$taxa[df$diffexpressed != "NO"]
#
# mycolors <- c("red", "red", "grey")
# names(mycolors) <- c("DOWN", "UP", "NO")
#
# ggplot(df, aes(x=summary.logFC, y=-log10(p.value), color=diffexpressed, label=delabel)) +
#   geom_point(alpha=0.6, size=2, show.legend = FALSE) +
#   geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype="longdash") + #
#   geom_hline(yintercept=-log10(0.05), col="black", linetype="longdash") +
#   geom_text_repel(color="black", size=3) +
#   scale_colour_manual(values = mycolors) +
#   theme_light() +
#   ggtitle("monocytes versus non-monocytes") +
#   ylab("-log10(P)") +
#   xlab("log2 fold change") +
#   ylim(0, 5)
#
#   # ylim(0, 12.5) +
#   # xlim(-2, 3) +
#
ggsave(snakemake@output[[1]], width=4, height=3, units = "in")
