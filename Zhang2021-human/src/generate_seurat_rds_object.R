library(Seurat)
library(stringr)
library(ggplot2)

# tumor.obj <- readRDS("output/tumor_cells.rds")

# cd45pos.read.matrix <- read.table("data/GSE160269_CD45pos_UMIs.txt", sep=" ", header = TRUE, row.names = 1)
cd45pos.read.matrix <- read.table(snakemake@input[[1]], sep=" ", header = TRUE, row.names = 1)
# cd45neg.read.matrix <- read.table("data/GSE160269_CD45neg_UMIs.txt", sep=" ", header = TRUE, row.names = 1)
cd45neg.read.matrix <- read.table(snakemake@input[[2]], sep=" ", header = TRUE, row.names = 1)
# meta.data <- read.table("data/metadata_for_host_transcriptome_analysis.tsv", sep="\t", header=TRUE, row.names = 1)
meta.data <- read.table(snakemake@input[[3]], sep="\t", header=TRUE, row.names = 1)

# use tumor.obj <- AddMetaData(tumor.obj, meta.data) to update metadata later
colnames(cd45pos.read.matrix) <- gsub("I", "CD45pos", colnames(cd45pos.read.matrix))
colnames(cd45pos.read.matrix) <- paste0(colnames(cd45pos.read.matrix), "-1")
colnames(cd45pos.read.matrix) <- gsub("-", ".", colnames(cd45pos.read.matrix))

colnames(cd45neg.read.matrix) <- gsub("E", "CD45neg", colnames(cd45neg.read.matrix))
colnames(cd45neg.read.matrix) <- paste0(colnames(cd45neg.read.matrix), "-1")
colnames(cd45neg.read.matrix) <- gsub("-", ".", colnames(cd45neg.read.matrix))

row.names(meta.data) <- gsub("-", ".", row.names(meta.data))

cd45pos.cell.names <- intersect(unique(colnames(cd45pos.read.matrix)), unique(rownames(meta.data)))

cd45pos.read.matrix <- cd45pos.read.matrix[, cd45pos.cell.names]
cd45pos.meta.data <- meta.data[cd45pos.cell.names,]

cd45pos.obj <- CreateSeuratObject(counts=cd45pos.read.matrix, meta.data=cd45pos.meta.data, project="cd45pos")

cd45neg.cell.names <- intersect(unique(colnames(cd45neg.read.matrix)), unique(rownames(meta.data)))

cd45neg.read.matrix <- cd45neg.read.matrix[, cd45neg.cell.names]
cd45neg.meta.data <- meta.data[cd45neg.cell.names,]

cd45neg.obj <- CreateSeuratObject(counts=cd45neg.read.matrix, meta.data=cd45neg.meta.data, project="cd45neg")

tumor.obj <- merge(cd45pos.obj, y=cd45neg.obj, add.cell.ids=c("CD45pos", "CD45neg"))
saveRDS(tumor.obj, "output/Zhang_raw_tumor.rds")
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

tumor.obj <- RunPCA(tumor.obj, features = VariableFeatures(object = tumor.obj))

# print(tumor.obj[["pca"]], dims = 1:5, nfeatures = 5)
tumor.obj <- JackStraw(tumor.obj, num.replicate = 100)
tumor.obj <- ScoreJackStraw(tumor.obj, dims = 1:20)
# JackStrawPlot(tumor.obj, dims = 1:15)

tumor.obj <- FindNeighbors(tumor.obj, dims = 1:15)
tumor.obj <- FindClusters(tumor.obj, resolution = 0.5)

saveRDS(tumor.obj, "output/Zhang_normalized_tumor.rds")

tumor.obj <- RunUMAP(tumor.obj, dims = 1:15)
p <- DimPlot(tumor.obj, reduction = "umap", group.by="celltype1")
ggsave("output/UMAP_plot.pdf", p)

saveRDS(tumor.obj, snakemake@output[[1]])
