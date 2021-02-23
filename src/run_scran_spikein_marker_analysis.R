library(scater)
library(scran)

source(snakemake@params[["spike_functions"]])

sce <- generate_sce(microbe.file=snakemake@input[[1]], spikein.file=snakemake@input[[2]],
  pdata.file=snakemake@input[[3]], celltype.col=snakemake@wildcards[["celltype"]])

celltype.col <- snakemake@wildcards[["celltype"]]
celltype.of.interest <- snakemake@wildcards[["celltype_of_interest"]]
celltype.comparison <- snakemake@wildcards[["celltype_comparison"]]

# now subset the cells by celltypes of interest if applicable
if (celltype.comparison != "all"){
  sce <- sce[, sce[[celltype.col]] %in% c(celltype.of.interest, celltype.comparison)]
}

min.proportion <- as.numeric(snakemake@wildcards[["minprop"]])
min.cpm <- 10

# get CPM using edgeR
SingleCellExperiment::cpm(sce) <- edgeR::cpm(counts(sce))

# if there is only one unique value in sce[[celltype.col]], this function returns a list
# if there is > 1 unique value in sce[[celltype.col]], this function returns a matrix
# get the proportion of cells with cpm > for all possible combinations of cell group and gene
propOverXbyGroup <- apply(
  cpm(sce),
  1,
  function(x){
  tapply(x, sce[[celltype.col]], function(x){sum(x > min.cpm) / length(x)})
})

if (length(unique(sce[[celltype.col]])) == 1){
  propOverXbyGroup <- t(as.matrix(propOverXbyGroup))
  rownames(propOverXbyGroup) <- unique(sce[[celltype.col]])
}

propOverXbyGroup <- reshape2::melt(
  propOverXbyGroup, varnames = c("group", "feature"), value.name = "Proportion")

# boolean for whether the proportion is > .25
propOverXbyGroup$PropOverXpct <- (propOverXbyGroup$Proportion > min.proportion)

# now we have a data.frame with four columns - group, feature, proportion and PropOver25pct

groupsOverXpct <- as.data.frame(with(
    propOverXbyGroup,
    tapply(PropOverXpct, feature, function(x){sum(x)}))
)
colnames(groupsOverXpct) <- "Groups"
groupsOverXpct$feature <- rownames(groupsOverXpct)

detectedFeatures <- names(which(tapply(
  propOverXbyGroup$Proportion, propOverXbyGroup$feature,
  function(x){sum(x > min.proportion) > 0})))

keepEndogenous <- (rownames(sce) %in% detectedFeatures)
sce <- sce[keepEndogenous]


celltype <- sce[[celltype.col]]

# filter so we are only comparing highly expressed OTUs to reduce FDR penalty
# agg <- aggregate(colnames(sce), by=list(sce[[celltype.col]]), FUN=length)
# min.cells <- min(agg$x)*.5
# print(min.cells)
# sce <- sce[rowSums(logcounts(sce) > 2) > min.cells, ]
# celltype.of.interest <- snakemake@wildcards[["celltype_of_interest"]]
groups <- sce[[celltype.col]]
lfc <- as.numeric(snakemake@wildcards[["lfc"]])
pval.type <- snakemake@wildcards[["pval.type"]]
block <- sce[[snakemake@wildcards[["block"]]]]
direction <- snakemake@wildcards[["direction"]]
# calculate markers using t-test
t.test.markers <- findMarkers(sce, groups=groups, lfc=lfc, pval.type=pval.type, block=block, direction=direction)
t.test.df <- t.test.markers[[celltype.of.interest]]
write.table(t.test.df, file=snakemake@output[[1]], sep="\t")
# calculate markers using wilcoxon rank sum test
wilcox.markers <- findMarkers(sce, groups=groups, test="wilcox", lfc=lfc, pval.type=pval.type, block=block, direction=direction)
wilcox.df <- wilcox.markers[[celltype.of.interest]]
write.table(wilcox.df, file=snakemake@output[[2]], sep="\t")
