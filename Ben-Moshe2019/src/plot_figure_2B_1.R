#library(EnhancedVolcano)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggrepel)
library(scales)


# input <- c("output/Pt0/GSM3454529/fisher-exact-Monocyte-Monocyte-nonMonocyte-phylum-PathSeq-Bacteria.tsv", "output/Pt0/GSM3454529/fisher-exact-Monocyte-Monocyte-nonMonocyte-class-PathSeq-Bacteria.tsv", "output/Pt0/GSM3454529/fisher-exact-Monocyte-Monocyte-nonMonocyte-order-PathSeq-Bacteria.tsv", "output/Pt0/GSM3454529/fisher-exact-Monocyte-Monocyte-nonMonocyte-family-PathSeq-Bacteria.tsv", "output/Pt0/GSM3454529/fisher-exact-Monocyte-Monocyte-nonMonocyte-genus-PathSeq-Bacteria.tsv", "output/Pt0/GSM3454529/fisher-exact-Monocyte-Monocyte-nonMonocyte-species-PathSeq-Bacteria.tsv")
my_data <- lapply(snakemake@input, read.table, header=TRUE)
#my_data <- lapply(input, read.table)
#print(my_data)

df <- do.call(rbind, my_data)
#print(df)

# taken from - https://stackoverflow.com/questions/47085514/simple-way-to-visualise-odds-ratios-in-r
# scales::pseudo_log_trans taken from https://stackoverflow.com/questions/40219639/how-to-deal-with-zero-in-log-plot

taxa <- rev(c("Proteobacteria", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae", "Salmonella", "Mycoplasma_wenyonii"))
df <- df[df$taxa %in% taxa,]
df$taxa <- factor(df$taxa, levels=taxa)
# levels(df$taxa) <- taxa

ggplot(df, aes(x = summary.odds.ratio, y = taxa)) +
  geom_vline(aes(xintercept = 1), size = .25, linetype="dashed") +
  geom_errorbarh(aes(xmax = or.ci.high, xmin=or.ci.low), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange") +
  theme_pubr(base_size=8) +
  ggtitle("monocytes versus non-monocytes (10x)") +
  ylab(NULL) +
  xlab("Odds ratio") +
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10))

ggsave(snakemake@output[[1]], width=3.5, height=3, units = "in")
