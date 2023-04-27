# library(scater)
# library(scran)
library(ggplot2)
library(ggpubr)
library(ggsci)


df <- read.table(snakemake@input[[1]], sep="\t", header=TRUE)
# df <- read.table("output/expected_actual_infected_cells_per_celltype1.tsv", sep="\t", header=TRUE)
n.rows <- dim(df)[1]
condition <- c(rep("expected", n.rows), rep("infected", n.rows))
celltype <- rep(df$celltype1, 2)
n.cells <- c(df$expected, df$infected)
df <- data.frame(condition, celltype, n.cells)


ggplot(df, aes(fill=condition, y=n.cells, x=celltype)) +
    geom_bar(position="dodge", stat="identity") +
    theme_pubr(base_size=6) + # , base_size= 10 to set font size
    ylab("number of cells") +
    theme(legend.title=element_blank()) +
    xlab(NULL) +
    scale_fill_manual(values=c("#4DBBD5FF", "#E64B35FF")) +
    ggtitle("All cells (Zhang2021)")
# scale_fill_manual() + # inspired by scale_color_npg - values=c("#00A087FF", "#4DBBD5FF","#E64B35FF")

ggsave(snakemake@output[[1]], width = 2.5, height = 2, units = "in")
