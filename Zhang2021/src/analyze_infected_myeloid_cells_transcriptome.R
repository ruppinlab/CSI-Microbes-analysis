library(Seurat)
library(stringr)
library(ggplot2)

# tumor.obj <- readRDS("output/tumor_cells.rds")

read.matrix <- read.table("data/GSE160269_CD45pos_UMIs.txt", sep=" ", header = TRUE, row.names = 1)
meta.data <- read.table("output/myeloid_cells_metadata.tsv", sep="\t", header=TRUE, row.names = 1)
# cd45pos_human_read_df.columns = cd45pos_human_read_df.columns.map(lambda x: x.replace("I", "CD45pos")+"-1")

# use tumor.obj <- AddMetaData(tumor.obj, meta.data) to update metadata later
colnames(read.matrix) <- gsub("I", "CD45pos", colnames(read.matrix))
colnames(read.matrix) <- paste0(colnames(read.matrix), "-1")
row.names(meta.data) <- gsub("-", ".", row.names(meta.data))
colnames(read.matrix) <- gsub("-", ".", colnames(read.matrix))

cell.names <- intersect(unique(colnames(read.matrix)), unique(rownames(meta.data)))

read.matrix <- read.matrix[, cell.names]
meta.data <- meta.data[cell.names,]

tumor.obj <- CreateSeuratObject(counts=read.matrix, meta.data=meta.data, project="tumor")

tumor.obj <- NormalizeData(tumor.obj, normalization.method = "LogNormalize", scale.factor = 10000)

tumor.obj <- FindVariableFeatures(tumor.obj, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(tumor.obj), 10)
plot1 <- VariableFeaturePlot(tumor.obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

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


Idents(object = tumor.obj) <- tumor.obj@meta.data$MyeloidSubtype
tam.obj <- subset(x = tumor.obj, idents = "TAM")
Idents(object = tam.obj) <- tam.obj@meta.data$infection
tam.markers <- FindMarkers(tam.obj, ident.1 = "infected", ident.2 = "uninfected", logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)

tam.markers <- FindMarkers(tam.obj, ident.1 = "infected", ident.2 = "uninfected")
head(tam.markers)
Idents(object = tam.obj) <- tam.obj@meta.data$microbe_level
tam.markers <- FindMarkers(tam.obj, ident.1 = "High", ident.2 = "None")
head(tam.markers)
tam.markers[c("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "APP"),]

tam.markers <- FindMarkers(tam.obj, ident.1 = "Low", ident.2 = "None")
head(tam.markers)
tam.markers[c("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "APP"),]

tam.markers <- FindMarkers(tam.obj, ident.1 = "High", ident.2 = "Low")
head(tam.markers)
tam.markers[c("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "APP"),]


Idents(object = tumor.obj) <- tumor.obj@meta.data$MyeloidSubtype
mono.obj <- subset(x = tumor.obj, idents = "Monocyte")
Idents(object = mono.obj) <- mono.obj@meta.data$infection
mono.markers <- FindMarkers(mono.obj, ident.1 = "infected", ident.2 = "uninfected")
head(mono.markers)
mono.markers[c("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "APP"),]

Idents(object = mono.obj) <- mono.obj@meta.data$microbe_level
mono.markers <- FindMarkers(mono.obj, ident.1 = "High", ident.2 = "None")
head(mono.markers)
mono.markers[c("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "APP"),]

Idents(object = mono.obj) <- mono.obj@meta.data$microbe_level
mono.markers <- FindMarkers(mono.obj, ident.1 = "Low", ident.2 = "None")
head(mono.markers)
mono.markers[c("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "APP"),]



markers <- FindMarkers(subset.tumor.obj, ident.1 = "infected", ident.2 = "uninfected", logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)



markers <- FindMarkers(tumor.obj, ident.1 = "infected", ident.2 = "uninfected") # ident.2 = "uninfected" # what about groupby="patient"?
markers <- FindMarkers(tumor.obj, ident.1 = "infection", groupby = "infection")

cluster2.markers <- FindMarkers(tumor.obj, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

VlnPlot(tumor.obj, features = c("HLA-A", "HLA-B", "HLA-C"))

# use this command to be DE genes for GSEA
infected.uninfected.markers <- FindMarkers(tumor.obj, ident.1 = "infected", ident.2 = "uninfected", logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)

GroupCorrelation(tumor.obj)
