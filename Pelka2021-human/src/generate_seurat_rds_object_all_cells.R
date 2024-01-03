library(Seurat)
library(stringr)
library(ggplot2)

# read.matrix <- Read10X_h5("data/GSE178341_crc10x_full_c295v4_submit.h5")
read.matrix <- Read10X_h5(snakemake@input[[1]])
# meta.data <- read.table("data/metadata_for_host_transcriptome_analysis.tsv", sep="\t", header=TRUE, row.names = 1)
meta.data <- read.table(snakemake@input[[2]], sep="\t", header=TRUE, row.names = 1)

cell.names <- intersect(unique(colnames(read.matrix)), unique(rownames(meta.data)))
read.matrix <- read.matrix[, cell.names]
meta.data <- meta.data[cell.names,]

tumor.obj <- CreateSeuratObject(counts=read.matrix, meta.data=meta.data, project="tumor")

# celltype.col <- snakemake@wildcards[["celltype_column"]]
# celltype <- snakemake@wildcards[["celltype"]]
# # subset to only desired cell-type
# Idents(object = tumor.obj) <- tumor.obj@meta.data[[celltype.col]]
# tumor.obj <- subset(tumor.obj, idents = celltype)

tumor.obj <- NormalizeData(tumor.obj, normalization.method = "LogNormalize", scale.factor = 10000)

tumor.obj <- FindVariableFeatures(tumor.obj, selection.method = "vst", nfeatures = 2000)

# top10 <- head(VariableFeatures(tumor.obj), 10)
# plot1 <- VariableFeaturePlot(tumor.obj)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2

all.genes <- rownames(tumor.obj)
tumor.obj <- ScaleData(tumor.obj, features = all.genes)

write.table(FetchData(object = tumor.obj, vars = c("CXCL8")), file="output/CXCL8_expression.tsv", sep="\t")

pcs <- c("PC_1", "PC_2", "PC_3", "PC_4", "PC_5", "PC_6", "PC_7", "PC_8", "PC_9", "PC_10", "PC_11", "PC_12", "PC_13", "PC_14", "PC_15")

write.table(FetchData(object = tumor.obj, vars = c()), file="output/CXCL8_expression.tsv", sep="\t")


FeaturePlot(tumor.obj, features = "CXCL8", cols = c("blue", "red"),
            split.by = c("celltype1", "infection"), 
            pt.size = 1, combine = FALSE)

tumor.obj <- RunPCA(tumor.obj, features = VariableFeatures(object = tumor.obj))

# print(tumor.obj[["pca"]], dims = 1:5, nfeatures = 5)
tumor.obj <- JackStraw(tumor.obj, num.replicate = 100)
tumor.obj <- ScoreJackStraw(tumor.obj, dims = 1:20)
# JackStrawPlot(tumor.obj, dims = 1:15)
tumor.obj <- FindNeighbors(tumor.obj, dims = 1:15)
tumor.obj <- FindClusters(tumor.obj, resolution = 0.5)

tumor.obj <- RunUMAP(tumor.obj, dims = 1:15)

p <- DimPlot(tumor.obj, reduction = "umap", group.by="celltype1")
ggsave("output/UMAP_plot.pdf", p)

saveRDS(tumor.obj, snakemake@output[[1]])
