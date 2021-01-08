library(Seurat)


pe.data <- Read10X(data.dir = snakemake@params[[1]])

se.data <- Read10X(data.dir = snakemake@params[[2]])

se.matrix <- as.data.frame(as.matrix(se.data))
pe.matrix <- as.data.frame(as.matrix(pe.data))


a <- setdiff(unique(colnames(pe.matrix)), unique(colnames(se.matrix)))
b <- setdiff(unique(colnames(se.matrix)), unique(colnames(pe.matrix)))
stopifnot(length(a) == 0)
stopifnot(length(b) == 0)
stopifnot(colnames(se.matrix) == colnames(pe.matrix))
m <- se.matrix + pe.matrix

write.table(m, snakemake@output[[1]], sep="\t")
