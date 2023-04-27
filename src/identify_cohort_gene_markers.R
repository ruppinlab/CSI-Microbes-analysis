library(Seurat)
library(stringr)
library(limma)


# tumor.obj <- readRDS("output/seurat_objects/celltype1_Myeloid.rds")
tumor.obj <- readRDS(snakemake@input[[1]])
# meta.data <- read.table("output/phenotypic_markers/infection_inclusive_2.tsv", sep="\t", header=TRUE, row.names = 1)
#meta.data <- read.table(snakemake@input[[2]], sep="\t", header=TRUE, row.names = 1)

#tumor.obj <- AddMetaData(tumor.obj, meta.data)

phenotype <- snakemake@wildcards[["phenotype"]]
ident1 <- snakemake@wildcards[["ident1"]]
ident2 <- snakemake@wildcards[["ident2"]]
Idents(object = tumor.obj) <- tumor.obj@meta.data[[phenotype]]
markers <- FindMarkers(tumor.obj, ident.1 = ident1, ident.2 = ident2)
write.table(markers, snakemake@output[[1]], sep="\t")
