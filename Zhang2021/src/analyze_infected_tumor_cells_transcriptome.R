library(Seurat)
library(stringr)
library(ggplot2)


read.matrix <- read.table("output/tumor_cells_human_microbe_reads.tsv", sep="\t", header = TRUE, row.names = 1)
# tumor.obj <- readRDS("output/tumor_cells.rds")

meta.data <- read.table("output/tumor_cells_metadata.tsv", sep="\t", header=TRUE, row.names = 1)
row.names(meta.data) <- gsub("-", ".", row.names(meta.data))
# use tumor.obj <- AddMetaData(tumor.obj, meta.data) to update metadata later

cell.names <- intersect(unique(colnames(read.matrix)), unique(rownames(meta.data)))
read.matrix <- read.matrix[, cell.names]
meta.data <- meta.data[cell.names,]

tumor.obj <- CreateSeuratObject(counts=read.matrix, meta.data=meta.data, project="tumor")

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
# tumor.obj <- JackStraw(tumor.obj, num.replicate = 100)
# tumor.obj <- ScoreJackStraw(tumor.obj, dims = 1:20)
# JackStrawPlot(tumor.obj, dims = 1:15)
tumor.obj <- FindNeighbors(tumor.obj, dims = 1:15)
tumor.obj <- FindClusters(tumor.obj, resolution = 0.5)

tumor.obj <- RunUMAP(tumor.obj, dims = 1:10)
DimPlot(tumor.obj, reduction = "umap")



Idents(object = tumor.obj) <- tumor.obj@meta.data$infection
markers <- FindMarkers(tumor.obj, ident.1 = "infected", ident.2 = "uninfected")
infected.uninfected.subset.markers <- FindMarkers(tumor.obj, ident.1 = "infected", ident.2 = "uninfected", logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)
subset.markers[c("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "APP"),]

# Pseudomonas
Idents(object = tumor.obj) <- tumor.obj@meta.data$Pseudomonas_infection
markers <- FindMarkers(tumor.obj, ident.1 = "infected", ident.2 = "uninfected")
infected.uninfected.subset.markers <- FindMarkers(tumor.obj, ident.1 = "infected", ident.2 = "uninfected", logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)
subset.markers[c("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "APP"),]

markers <- FindMarkers(tumor.obj, ident.1 = "infected", ident.2 = "uninfected") # ident.2 = "uninfected" # what about groupby="patient"?
markers <- FindMarkers(tumor.obj, ident.1 = "infection", groupby = "infection")

cluster2.markers <- FindMarkers(tumor.obj, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

VlnPlot(tumor.obj, features = c("HLA-A", "HLA-B", "HLA-C"))

# use this command to be DE genes for GSEA
infected.uninfected.markers <- FindMarkers(tumor.obj, ident.1 = "infected", ident.2 = "uninfected", logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)



df <- FetchData(object = tumor.obj, vars = c("infection", "Pseudomonas_infection", "Sphingomonas_infection", "Mycoplasma_infection", "HLA-A", "HLA-B", "HLA-C", "B2M", "HLA-DQA1", "HLA-DQA2", "HLA-DOA", "HLA-DPB1", "CD74", "HLA-DRA", "HLA-DPA1", "HLA-DRB5", "HLA-DRB1", "HLA-DMA", "HLA-DMB", "HLA-DQB1", "HLA-DOB", "HLA-DQB2"))
colnames(df) <- gsub("-", ".", colnames(df))
df$HLA.I <- df$HLA.A + df$HLA.B + df$HLA.C + df$B2M
df$HLA.II <- df$HLA.DQA1 + df$HLA.DQA2 + df$HLA.DOA + df$HLA.DPB1 + df$CD74 + df$HLA.DRA + df$HLA.DPA1 + df$HLA.DRB5 + df$HLA.DRB1 + df$HLA.DMA + df$HLA.DMB + df$HLA.DQB1 + df$HLA.DOB + df$HLA.DQB2
