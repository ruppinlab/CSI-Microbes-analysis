library(EnhancedVolcano)

# markers <- read.table("output/DE_genes/celltype1_Myeloid_infection_inclusive_2.tsv", sep="\t")
markers <- read.table(snakemake@input[[1]], sep="\t")

EnhancedVolcano(markers,
  lab = rownames(markers),
  x = 'avg_log2FC',
  y = 'p_val',
  xlim = c(-2, 2),
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  FCcutoff = .25,
  labSize=1.5,
  pointSize=.5,
  ylab=NULL,
  xlab=NULL,
  legendLabSize=NULL,
  captionLabSize=NULL,
  axisLabSize=6,
  legendLabels=NULL,
  legendPosition="none",
  legendIconSize=NULL,
  pCutoff = 10e-6)
#   selectLab=c("HLA-DPA1", "HLA-DOA", "HLA-DPB1", "HLA-DQA1", "HLA-DOB", "HLA-DQB1", "HLA-DRB1"),
#   drawConnectors=TRUE,
#   max.overlaps=50,


# ggsave("output/plots/celltype1_Epi_infection_inclusive_2.pdf")
ggsave(snakemake@output[[1]], units="in", width=5, height=4)


EnhancedVolcano(markers,
  lab = rownames(markers),
  x = 'avg_log2FC',
  y = 'p_val',
  xlim = c(-2.5, 2.5),
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  FCcutoff = .5,
  labSize=1,
  pointSize=.5,
  ylab=NULL,
  xlab=NULL,
  legendLabSize=NULL,
  captionLabSize=NULL,
  axisLabSize=6,
  legendLabels=NULL,
  legendPosition="none",
  legendIconSize=NULL,
  pCutoff = 10e-6)
