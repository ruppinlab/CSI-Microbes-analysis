library(scater)
library(scuttle)
library(ggplot2)


sce <- readRDS("output/sce_filtered_obj.rds")
sce$celltype <- factor(colData(sce)$Cell_Type)
sce$condition <- factor(colData(sce)$Condition)
# read in file with Fusobacterium abundance
meta <- read.table("output/meta_data_for_UMAP.tsv", sep="\t", header=TRUE, row.names=1)
Fusobacterium <- meta[colnames(sce), "Fusobacterium"]
colData(sce) <- cbind(colData(sce), Fusobacterium)

fig <- cowplot::plot_grid(ncol = 4,
  plotUMAP(sce, colour_by="label"),
  plotUMAP(sce, colour_by="celltype"),
  plotUMAP(sce, colour_by="condition"),
  plotUMAP(sce, colour_by="Fusobacterium")
  )

ggsave("output/plots/UMAP_plot.svg", fig, width=8, height=2)
