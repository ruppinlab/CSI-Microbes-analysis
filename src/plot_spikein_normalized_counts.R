library(scater)
library(scran)

source(snakemake@params[["spike_functions"]])

sce <- generate_sce(microbe.file=snakemake@input[[1]], spikein.file=snakemake@input[[2]],
  pdata.file=snakemake@input[[3]], celltype.col=snakemake@wildcards[["celltype"]])

tax.map <- read.table(snakemake@input[[4]], sep="\t", header=TRUE)
species.tax.map <- tax.map[tax.map$taxa_level == snakemake@wildcards[["tax_level"]],]

# rename rownames from tax_id to names
rownames(sce) <- lapply(rownames(sce), function(x) tax.map[tax.map$tax_id == x, "name"])

p <- plotExpression(sce, x=snakemake@wildcards[["celltype"]], features=snakemake@wildcards[["microbe"]]) +
    theme(axis.text.x = element_text(angle = 90))
p <- p + ylab("log2(spike-in normalized abundance + 1)")
p <- p + ggtitle(snakemake@wildcards[["microbe"]])

ggsave(snakemake@output[[1]], p)
