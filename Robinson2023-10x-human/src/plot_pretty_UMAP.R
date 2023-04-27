library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
library(ggplot2)
# install.packages("scCustomize")
library(scCustomize)

set.seed(0)

exposure.annotate <- function(orig.ident) {
  if (orig.ident == "uninfected"){
    return ("Unexposed")
  }
  if (orig.ident == "HK"){
    return ("Heat-killed")
  }
  if (orig.ident == "live") {
    return ("Live")
  }
}




combined <- readRDS("output/combined_seurat_obj.rds")
combined@meta.data$exposure <- sapply(combined@meta.data$orig.ident, exposure.annotate)
combined@meta.data$Fuso_log <- log2(combined@meta.data$Fusobacterium+1)
combined <- FindNeighbors(combined, dims = 1:15)
combined <- FindClusters(combined, resolution = .5)
combined <- RunUMAP(combined, dims = 1:15)
DimPlot(combined, reduction = "umap")

# let's plot EPCAM, CD3E and something and label HCT116, Jurkat and THP1
# epcam <- WhichCells(object = combined, expression = EPCAM > 1.5)
# cd3e <- WhichCells(combined, expression = CD3E > 1.5)
# cd64 <- WhichCells(object = combined, expression = FCGR1A > .75)
# cells <- list(EPCAM = epcam, CD3E = cd3e, CD64 = cd64)
# Cell_Highlight_Plot(combined, cells_highlight=cells, highlight_color=JCO_Four()[1:3])
# Cell_Highlight_Plot(combined, cells_highlight=cells, highlight_color=c("#0073C2FF", "#EFC000FF", "#CD534CFF"))
# Fusobacterium plot - would be great to better show lower difference
pal <- viridis(n = 10, option = "C", direction = -1)

fig <- cowplot::plot_grid(ncol = 4,
  DimPlot(combined, group.by = "celltype", shuffle=TRUE) + labs(title=NULL, subtitle=NULL),
  DimPlot(combined) + labs(title=NULL, subtitle=NULL),
  DimPlot(combined, group.by = "exposure", shuffle=TRUE) + labs(title=NULL, subtitle=NULL),
  FeaturePlot_scCustom(combined, features = "Fuso_log", colors_use=pal) + labs(title=NULL, subtitle=NULL)
  )
ggsave("output/plots/UMAP_all_cells_with_legend.svg", scale=2, height=2, width=8.5)

fig <- cowplot::plot_grid(ncol = 4,
  DimPlot(combined, group.by = "celltype", shuffle=TRUE) + labs(title=NULL, subtitle=NULL) + theme(legend.position="none"),
  DimPlot(combined) + labs(title=NULL, subtitle=NULL) + theme(legend.position="none"),
  DimPlot(combined, group.by = "exposure", shuffle=TRUE) + labs(title=NULL, subtitle=NULL) + theme(legend.position="none"),
  FeaturePlot_scCustom(combined, features = "Fuso_log", colors_use=pal) + labs(title=NULL, subtitle=NULL) + theme(legend.position="none")
  )
ggsave("output/plots/UMAP_all_cells_no_legend.svg", scale=2, height=2, width=6)


fig <- cowplot::plot_grid(ncol = 2,
  DimPlot(combined, group.by = "celltype1", shuffle=TRUE) + labs(title=NULL, subtitle=NULL),
  FeaturePlot_scCustom(combined, features = "log_microbial_UMIs", colors_use=pal) + labs(title=NULL, subtitle=NULL)
  )
ggsave("Pelka2021/output/plots/UMAP_all_cells_with_legend.svg", scale=2, height=2, width=4)

fig <- cowplot::plot_grid(ncol = 2,
  DimPlot(combined, group.by = "celltype1", shuffle=TRUE) + labs(title=NULL, subtitle=NULL) + theme(legend.position="none"),
  FeaturePlot_scCustom(combined, features = "log_microbial_UMIs", colors_use=pal) + labs(title=NULL, subtitle=NULL) + theme(legend.position="none")
  )
ggsave("Pelka2021/output/plots/UMAP_all_cells_no_legend.svg", scale=2, height=2, width=4)


# fig <- cowplot::plot_grid(ncol = 4,
#   DimPlot(combined, group.by = "celltype", shuffle=TRUE) + labs(title=NULL, subtitle=NULL) + theme(legend.position="none"),
#   DimPlot(combined) + labs(title=NULL, subtitle=NULL) + theme(legend.position="none"),
#   DimPlot(combined, group.by = "exposure", shuffle=TRUE) + labs(title=NULL, subtitle=NULL) + theme(legend.position="none"),
#   FeaturePlot_scCustom(combined, features = "Fuso_log", colors_use=pal) + labs(title=NULL, subtitle=NULL) + theme(legend.position="none")
#   )
# ggsave("output/plots/UMAP_all_cells_no_legend.svg", scale=2, height=2, width=6)
