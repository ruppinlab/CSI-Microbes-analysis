library(dplyr)
library(Seurat)
library(patchwork)
suppressMessages(require(DoubletFinder))

set.seed(0)


data <- Read10X(data.dir = "raw/P1-SCAF2961_1_Uninfected/outs/filtered_feature_bc_matrix")
uninfected <- CreateSeuratObject(counts = data, project="SCAF2961_1_Uninfected", min.cells = 3, min.features=200)
data <- Read10X(data.dir = "raw/P1-SCAF2962_2_HK/outs/filtered_feature_bc_matrix")
HK <- CreateSeuratObject(counts = data, project="SCAF2962_2_HK", min.cells = 3, min.features=200)
data <- Read10X(data.dir = "raw/P1-SCAF2963_3_Live/outs/filtered_feature_bc_matrix")
live <- CreateSeuratObject(counts = data, project="SCAF2963_3_Live", min.cells = 3, min.features=200)

combined <- merge(uninfected, y=c(HK, live), add.cell.ids=c("SCAF2961_1_Uninfected", "SCAF2962_2_HK", "SCAF2963_3_Live"))

combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")

# VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

combined <- subset(combined, subset = percent.mt < 15 & nFeature_RNA > 3000)

combined <- NormalizeData(combined, normalization.method="LogNormalize", scale.factor = 10000)

combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(combined), 10)
plot1 <- VariableFeaturePlot(combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2
# scale the data
all.genes <- rownames(combined)
combined <- ScaleData(combined, features = all.genes)

# perform linear dimensional reduction
combined <- RunPCA(combined, features = VariableFeatures(object = combined))
# Examine and visualize PCA results a few different ways
#print(combined[["pca"]], dims = 1:5, nfeatures = 5)
#VizDimLoadings(combined, dims = 1:2, reduction = "pca")
#DimPlot(combined, reduction = "pca")
combined <- JackStraw(combined, num.replicate = 100)
combined <- ScoreJackStraw(combined, dims = 1:20)

#JackStrawPlot(combined, dims = 1:15)
#ElbowPlot(combined)
# cluster the cells
combined <- FindNeighbors(combined, dims = 1:15)
combined <- FindClusters(combined, resolution = 0.1)

combined <- RunUMAP(combined, dims = 1:15)
DimPlot(combined, reduction = "umap")
DimPlot(combined, reduction = "umap", group.by = "orig.ident")

# predict doublets
nExp <- round(ncol(combined) * 0.04)
combined <- doubletFinder_v3(combined, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

# visualize the doublets
DF.name = colnames(combined@meta.data)[grepl("DF.classification", colnames(combined@meta.data))]
cowplot::plot_grid(ncol = 2, DimPlot(combined, group.by = "orig.ident") + NoAxes(),
    DimPlot(combined, group.by = DF.name) + NoAxes())
VlnPlot(combined, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
# remove the doublets
combined <- combined[, combined@meta.data[, DF.name] == "Singlet"]
# re-run UMAP
combined <- FindNeighbors(combined, dims = 1:15)
combined <- FindClusters(combined, resolution = 0.1)
combined <- RunUMAP(combined, dims = 1:15)
DimPlot(combined, reduction = "umap")
DimPlot(combined, reduction = "umap", group.by = "orig.ident")

# cluster 4 looks like doublets of THP1 and Jurkat
VlnPlot(combined, features = "nFeature_RNA")
# remove cluster 4, which consists of only 39 cells
combined <- subset(combined, idents = 4, invert=TRUE)
# let's rename orig.idents to condition
combined@meta.data$condition <- combined@meta.data$orig.ident
# let's add cell-type annotations
celltype.annotate <- function(cluster) {
  if (cluster == 0){
    return ("HCT116")
  }
  if (cluster == 1){
    return ("Jurkat")
  }
  if (cluster %in% c(2,3)) {
    return ("THP1")
  }
}

combined@meta.data$celltype <- sapply(Idents(combined), celltype.annotate)

exposure.annotate <- function(orig.ident) {
  if (orig.ident == "SCAF2961_1_Uninfected"){
    return ("Unexposed")
  }
  if (orig.ident == "SCAF2962_2_HK"){
    return ("Heat-killed")
  }
  if (orig.ident == "SCAF2963_3_Live") {
    return ("Live")
  }
}

combined@meta.data$exposure <- sapply(combined@meta.data$orig.ident, exposure.annotate)


saveRDS(combined, "output/combined_3p_cells.rds")

# cowplot::plot_grid(ncol = 3,
#   DimPlot(combined,) + NoAxes(),
#   DimPlot(combined, group.by = "condition") + NoAxes(),
#   DimPlot(combined, group.by = "celltype") + NoAxes())
# ggsave("output/plots/UMAP_all_cells.pdf", height=7, width=18)
#
# # let's perform UMAP on specific cell-types
# thp1 <- subset(combined, idents=c(2,3))
# thp1 <- FindNeighbors(thp1, dims = 1:15)
# thp1 <- FindClusters(thp1, resolution = 0.1)
# thp1 <- RunUMAP(thp1, dims = 1:15)
# cowplot::plot_grid(ncol = 2, DimPlot(thp1) + NoAxes(),
#     DimPlot(thp1, group.by = "condition") + NoAxes())
# ggsave("output/plots/UMAP_THP1_cells.pdf", height=7, width=12)
#
# hct116 <- subset(combined, idents=c(0))
# hct116 <- FindNeighbors(hct116, dims = 1:15)
# hct116 <- FindClusters(hct116, resolution = 0.1)
# hct116 <- RunUMAP(hct116, dims = 1:15)
# cowplot::plot_grid(ncol = 2, DimPlot(hct116) + NoAxes(),
#     DimPlot(hct116, group.by = "condition") + NoAxes())
# ggsave("output/plots/UMAP_HCT116_cells.pdf", height=7, width=12)
#
# jurkat <- subset(combined, idents=c(1))
# jurkat <- FindNeighbors(jurkat, dims = 1:15)
# jurkat <- FindClusters(jurkat, resolution = 0.1)
# jurkat <- RunUMAP(jurkat, dims = 1:15)
# cowplot::plot_grid(ncol = 2, DimPlot(jurkat) + NoAxes(),
#     DimPlot(jurkat, group.by = "condition") + NoAxes())
# ggsave("output/plots/UMAP_Jurkat_cells.pdf", height=7, width=12)

write.table(combined@meta.data, file="output/filtered-3p-units.tsv", sep="\t")
