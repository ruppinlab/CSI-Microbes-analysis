library(scater)
library(scran)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(stringr)

sce <- readRDS(snakemake@input[[1]])

gene.of.interest <- snakemake@wildcards[["gene_of_interest"]]

sce$patient <- factor(sce$patient, levels=c("TH067", "TH103", "TH171", "TH179", "TH220", "TH225", "TH226", "TH227", "TH231", "TH236", "TH238", "TH248", "TH266"))

groups <- ifelse(sce$patient %in% c("TH231", "TH236", "TH238", "TH266"), "infected", "uninfected")

df <- data.frame(logcounts=logcounts(sce)[gene.of.interest, ], patient=sce$patient, group=groups)

# reference - https://stackoverflow.com/questions/4973898/combining-paste-and-expression-functions-in-plot-labels
ylabel <- bquote("Expression " * log[2] * "(spike-normalized " * .(gene.of.interest) * ")")

fig.b <- ggboxplot(df, x = "patient", y = "logcounts", color = "group", palette = c("#00AFBB", "#E7B800"), add = "jitter") +
rotate_x_text() +
ylab(ylabel)

fig.a <- ggboxplot(df, x = "group", y = "logcounts", color = "group", palette = c("#00AFBB", "#E7B800"), add = "jitter") +
ylab(ylabel)

ggarrange(fig.a, fig.b, ncol=2, nrow=1, labels=c("A", "B"))


ggsave(snakemake@output[[1]], width=16, height=6)
