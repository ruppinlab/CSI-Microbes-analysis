library(dplyr)
library(Seurat)
library(patchwork)
suppressMessages(require(DoubletFinder))


data <- Read10X_h5("raw/P1-SCAF2965_5_Live/outs/raw_feature_bc_matrix.h5")
live <- CreateSeuratObject(counts = data, project="live", min.cells = 3, min.features=200)

live[["percent.mt"]] <- PercentageFeatureSet(live, pattern = "^MT-")

VlnPlot(live, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# weird distribution for nFeature_RNA
live <- subset(live, subset = percent.mt < 10 & nFeature_RNA > 4000) # 4000 may be too high

live <- NormalizeData(live, normalization.method="LogNormalize", scale.factor = 10000)

live <- FindVariableFeatures(live, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(live), 10)
plot1 <- VariableFeaturePlot(live)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2
# scale the data
all.genes <- rownames(live)
live <- ScaleData(live, features = all.genes)

# perform linear dimensional reduction
live <- RunPCA(live, features = VariableFeatures(object = live))
# Examine and visualize PCA results a few different ways
# print(live[["pca"]], dims = 1:5, nfeatures = 5)
# VizDimLoadings(live, dims = 1:2, reduction = "pca")
# DimPlot(live, reduction = "pca")
live <- JackStraw(live, num.replicate = 100)
live <- ScoreJackStraw(live, dims = 1:20)

JackStrawPlot(live, dims = 1:15)
ElbowPlot(live)
# cluster the cells
live <- FindNeighbors(live, dims = 1:15)
live <- FindClusters(live, resolution = 0.1)

live <- RunUMAP(live, dims = 1:15)
DimPlot(live, reduction = "umap")


nExp <- round(ncol(live) * 0.04)
live <- doubletFinder_v3(live, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

DF.name = colnames(live@meta.data)[grepl("DF.classification", colnames(live@meta.data))]
# cowplot::plot_grid(ncol = 2, DimPlot(live, group.by = "orig.ident") + NoAxes(),
#     DimPlot(live, group.by = DF.name) + NoAxes())

VlnPlot(live, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
live = live[, live@meta.data[, DF.name] == "Singlet"]
live <- FindNeighbors(live, dims = 1:15)
live <- FindClusters(live, resolution = 0.1)

live <- RunUMAP(live, dims = 1:15)
DimPlot(live, reduction = "umap")

celltype.annotate <- function(cluster) {
  if (cluster %in% c(0, 3, 4)){
    return ("HCT116")
  }
  if (cluster == 2){
    return ("Jurkat")
  }
  if (cluster == 1) {
    return ("THP1")
  }
}

live@meta.data$celltype <- sapply(Idents(live), celltype.annotate)


write.table(live@meta.data, file="output/live5_units.tsv", sep="\t")
