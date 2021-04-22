library(dplyr)
library(scater)
library(scran)
#library(ggplot2)
#library(EnhancedVolcano)

counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
# counts <- read.table("output/Pt0_genus_PathSeq_microbe_reads.tsv", sep="\t", header=TRUE, row.names=1)
# output/61_species_PathSeq_metadata.tsv
pdata <- read.table(snakemake@input[[2]], sep="\t", header = TRUE, row.names=1)
row.names(pdata) <- gsub("-", ".", row.names(pdata))
# pdata <- read.table("output/Pt0_genus_PathSeq_metadata.tsv", sep="\t", header = TRUE, row.names=1)
tax.map <- read.table(snakemake@input[[3]], sep="\t", header=TRUE)

#
celltype.col <- snakemake@wildcards[["celltype"]]
celltype.of.interest <- snakemake@wildcards[["celltype_of_interest"]]
celltype.comparison <- snakemake@wildcards[["celltype_comparison"]]
print(dim(counts))
print(class(counts))
# if there are no reads, then write empty data.frame
if (nrow(counts) == 0){
  print("There are no reads at all. Writing empty data.frame")
  df <- data.frame(Top=character(), p.value=double(), FDR=double(), summary.logFC=double())
  if (celltype.comparison != "all"){
    n <- paste0("LogFC.", celltype.comparison)
    df[n] = double()
  } else {
    cols <- setdiff(unique(pdata[[celltype.col]]), celltype.of.interest)
    for (c in cols)
    {
      n <- paste0("LogFC.", c)
      df[n] <- double()
    }
  }
  write.table(df, file=snakemake@output[[1]], sep="\t")
  quit()
}

sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData=pdata)

# now subset the cells by celltypes of interest if applicable
if (celltype.comparison != "all"){
  sce <- sce[, sce[[celltype.col]] %in% c(celltype.of.interest, celltype.comparison)]
}

# set the minimum proportion of cells of any cell-type that must have at least one UMI
min.proportion <- as.numeric(snakemake@wildcards[["minprop"]])

# get the proportion of cells with count > 0 for all possible combinations of cell group and gene
propOver0byGroup <- apply(
  counts(sce),
  1,
  function(x){
  tapply(x, sce[[celltype.col]], function(x){sum(x > 0) / length(x)})
})

if (length(unique(sce[[celltype.col]])) == 1){
  propOver0byGroup <- t(as.matrix(propOver0byGroup))
  rownames(propOver0byGroup) <- unique(sce[[celltype.col]])
}

propOver0byGroup <- reshape2::melt(
  propOver0byGroup, varnames = c("group", "feature"), value.name = "Proportion")

# boolean for whether the proportion is > 2%
propOver0byGroup$PropOverXpct <- (propOver0byGroup$Proportion > min.proportion)

# now we have a data.frame with four columns - group, feature, proportion and PropOver25pct

groupsOverXpct <- as.data.frame(with(
    propOver0byGroup,
    tapply(PropOverXpct, feature, function(x){sum(x)}))
)

colnames(groupsOverXpct) <- "Groups"
groupsOverXpct$feature <- rownames(groupsOverXpct)

detectedFeatures <- names(which(tapply(
  propOver0byGroup$Proportion, propOver0byGroup$feature,
  function(x){sum(x > min.proportion) > 0})))

keepEndogenous <- (rownames(sce) %in% detectedFeatures)
sce <- sce[keepEndogenous]


celltype <- sce[[celltype.col]]

groups <- celltype
lfc <- as.numeric(snakemake@wildcards[["lfc"]])
pval.type <- snakemake@wildcards[["pvaltype"]]
block <- sce[[snakemake@wildcards[["block"]]]]
direction <- snakemake@wildcards[["direction"]]

binom <- findMarkers(sce, test="binom", groups=groups, pval.type=pval.type, block=block, assay.type="counts", direction=direction)

df <- binom[[celltype.of.interest]]
df$taxa <- lapply(rownames(df), function(x) tax.map[tax.map$tax_id == x, "name"])
write.table(df, file=snakemake@output[[1]], sep="\t")

#EnhancedVolcano(df, lab=df$taxa, x = "summary.logFC", y = "p.value", pCutoff = .01)
#ggsave(snakemake@output[[2]])
