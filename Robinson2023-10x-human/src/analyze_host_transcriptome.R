library(Seurat)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggsci)


# tumor.obj <- readRDS("output/tumor_cells.rds")

uninf.read.matrix <- Read10X_h5("data/SCAF2961_1_Uninfected_filtered_feature_bc_matrix.h5")
colnames(uninf.read.matrix) <- lapply(colnames(uninf.read.matrix), function(x){paste0("SCAF2961_1_Uninfected-", x)})

HK.read.matrix <- Read10X_h5("data/SCAF2962_2_HK_filtered_feature_bc_matrix.h5")
colnames(HK.read.matrix) <- lapply(colnames(HK.read.matrix), function(x){paste0("SCAF2962_2_HK-", x)})

live.read.matrix <- Read10X_h5("data/SCAF2963_3_Live_filtered_feature_bc_matrix.h5")
colnames(live.read.matrix) <- lapply(colnames(live.read.matrix), function(x){paste0("SCAF2963_3_Live-", x)})

CTP.read.matrix <- Read10X_h5("data/SCAF2964_4_HCT116_Cell_Trace_pos_filtered_feature_bc_matrix.h5")
colnames(CTP.read.matrix) <- lapply(colnames(CTP.read.matrix), function(x){paste0("SCAF2964_4_HCT116_Cell_Trace_pos-", x)})

read.matrix <- cbind(uninf.read.matrix, HK.read.matrix, live.read.matrix, CTP.read.matrix)

meta.data <- read.table("output/HCT116_DE_metadata.tsv", sep="\t", header=TRUE, row.names = 1)

cell.names <- intersect(unique(row.names(meta.data)), unique(colnames(read.matrix)))
meta.data <- meta.data[cell.names,]
read.matrix <- read.matrix[, cell.names]
# use tumor.obj <- AddMetaData(tumor.obj, meta.data) to update metadata later
# row.names(meta.data) <- gsub("-", ".", row.names(meta.data))

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

print(tumor.obj[["pca"]], dims = 1:5, nfeatures = 5)
tumor.obj <- JackStraw(tumor.obj, num.replicate = 100)
tumor.obj <- ScoreJackStraw(tumor.obj, dims = 1:20)
JackStrawPlot(tumor.obj, dims = 1:15)
tumor.obj <- FindNeighbors(tumor.obj, dims = 1:15)
tumor.obj <- FindClusters(tumor.obj, resolution = 0.5)

tumor.obj <- RunUMAP(tumor.obj, dims = 1:10)
DimPlot(tumor.obj, reduction = "umap")

df <- FetchData(object = tumor.obj, vars = c("exposure", "HLA-A", "HLA-B", "HLA-C", "B2M", "HLA-DQA1", "HLA-DQA2", "HLA-DOA", "HLA-DPB1", "CD74", "HLA-DRA", "HLA-DPA1", "HLA-DRB5", "HLA-DRB1", "HLA-DMA", "HLA-DMB", "HLA-DQB1", "HLA-DOB", "HLA-DQB2"))
colnames(df) <- gsub("-", ".", colnames(df))
df$HLA.I <- df$HLA.A + df$HLA.B + df$HLA.C + df$B2M
df$HLA.II <- df$HLA.DQA1 + df$HLA.DQA2 + df$HLA.DOA + df$HLA.DPB1 + df$CD74 + df$HLA.DRA + df$HLA.DPA1 + df$HLA.DRB5 + df$HLA.DRB1 + df$HLA.DMA + df$HLA.DMB + df$HLA.DQB1 + df$HLA.DOB + df$HLA.DQB2

df$exposure <- factor(df$exposure, levels=c("unexposed", "HK", "live"))
levels(df$exposure) <- c("Unexposed", "Heat-killed", "Live")

fig.hla.a <- ggplot(df, alpha=.25, aes(x=exposure, fill=exposure, y=HLA.A)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position=position_jitterdodge(), alpha=0.1, size=2, aes(color=exposure)) +
  theme_pubr(base_size=11) + # , base_size= 10 to set font size
  ylab("HLA-A expression") + xlab(NULL) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())# + scale_color_npg()
  # scale_color_manual(values=c("#00A087FF", "#E64B35FF")) + # inspired by scale_color_npg
  # ggtitle("Salmonella abundance")

fig.hla.b <- ggplot(df, alpha=.25, aes(x=exposure, fill=exposure, y=HLA.B)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position=position_jitterdodge(), alpha=0.1, size=2, aes(color=exposure)) +
  theme_pubr(base_size=11) + # , base_size= 10 to set font size
  ylab("HLA-B expression") + xlab(NULL) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

fig.hla.c <- ggplot(df, alpha=.25, aes(x=exposure, fill=exposure, y=HLA.C)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position=position_jitterdodge(), alpha=0.1, size=2, aes(color=exposure)) +
  theme_pubr(base_size=11) + # , base_size= 10 to set font size
  ylab("HLA-C expression") + xlab(NULL) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

fig.b2m <- ggplot(df, alpha=.25, aes(x=exposure, fill=exposure, y=B2M)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position=position_jitterdodge(), alpha=0.1, size=2, aes(color=exposure)) +
  theme_pubr(base_size=11) + # , base_size= 10 to set font size
  ylab("B2M expression") + xlab(NULL) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

fig.hla.II <- ggplot(df, alpha=.25, aes(x=exposure, fill=exposure, y=HLA.II)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position=position_jitterdodge(), alpha=0.1, size=2, aes(color=exposure)) +
  theme_pubr(base_size=11) + # , base_size= 10 to set font size
  ylab("HLA-II expression") + xlab(NULL) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

Idents(object = tumor.obj) <- tumor.obj@meta.data$exposure
live.unexposed.markers <- FindMarkers(tumor.obj, ident.1 = "live", ident.2 = "unexposed", logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)

# exposed.HK.markers <- FindMarkers(tumor.obj, ident.1 = "live", ident.2 = "HK")
live.HK.markers <- FindMarkers(tumor.obj, ident.1 = "live", ident.2 = "HK", logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)

live.exposed.tumor.obj <- subset(tumor.obj, idents = "live")
Idents(object = live.exposed.tumor.obj) <- live.exposed.tumor.obj@meta.data$infection
# subset.markers <- FindMarkers(subset.tumor.obj, ident.1 = "infected", ident.2 = "uninfected", logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)
live.exposed.infected.uninfected.markers <- FindMarkers(live.exposed.tumor.obj, ident.1 = "infected", ident.2 = "uninfected")
live.exposed.infected.uninfected.markers[c("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "APP"),]



markers <- FindMarkers(tumor.obj, ident.1 = "infected", ident.2 = "uninfected") # ident.2 = "uninfected" # what about groupby="patient"?
markers <- FindMarkers(tumor.obj, ident.1 = "infection", groupby = "infection")

VlnPlot(tumor.obj, features = c("HLA-A", "HLA-B", "HLA-C"))

# use this command to be DE genes for GSEA
infected.uninfected.markers <- FindMarkers(tumor.obj, ident.1 = "infected", ident.2 = "uninfected", logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)
