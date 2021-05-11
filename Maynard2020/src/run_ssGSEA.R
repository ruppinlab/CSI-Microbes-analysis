library(GSVA)

df <- read.table(snakemake@input[[1]], header=TRUE)
tcga.df <- read.table(snakemake@input[[2]], header=TRUE, row.names=1)

infected.down.gene.set <- df[df["FDR"] < .05 & df["summary.AUC"] < .5,]$gene
infected.up.gene.set <- df[df["FDR"] < .05 & df["summary.AUC"] > .5,]$gene
gene.set.list <- list(infected.dn=infected.down.gene.set, infected.up=infected.up.gene.set)

out <- gsva(as.matrix(tcga.df), gene.set.list, method="ssgsea")

write.table(out, snakemake@output[[1]])
