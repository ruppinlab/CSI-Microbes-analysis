library(scater)
library(scran)


human.counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
# human.counts <- read.table("output/star/Pt0/human_genes.tsv", sep="\t", header=TRUE, row.names=1)
spikein.counts <- read.table(snakemake@input[[2]], sep="\t", header=TRUE, row.names=1)
# spikein.counts <- read.table("output/star/Pt0/spike_ins.tsv", sep="\t", header=TRUE, row.names=1)
column.names <- intersect(unique(colnames(human.counts)), unique(colnames(spikein.counts)))
# remove samples with no spike-ins
column.names <- column.names[!(column.names %in% c("B1_B003648", "P5_B000420"))]
spikein.counts <- spikein.counts[, column.names]
human.counts <- human.counts[, column.names]
# row.names=2 because 2nd row is cell
pdata <- read.table(snakemake@input[[3]], sep="\t", header = TRUE, row.names=2)
#pdata <- read.table("data/units.tsv", sep="\t", header = TRUE, row.names=2)
pdata <- pdata[column.names, ]
row.names(pdata) <- gsub("-", ".", row.names(pdata))

# create SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(counts = as.matrix(human.counts)), colData=pdata)
spike_se <- SummarizedExperiment(list(counts=spikein.counts))
altExp(sce, "spike") <- spike_se
sce <- computeSpikeFactors(sce, "spike")
sce <- logNormCounts(sce)

celltype.col <- snakemake@wildcards[["celltype"]]
#celltype.col <- "infected"
celltype.of.interest <- snakemake@wildcards[["celltype_of_interest"]]
#celltype.of.interest <- "infected"
celltype.comparison <- snakemake@wildcards[["celltype_comparison"]]
#celltype.comparison <- "bystander"

# now subset the cells by celltypes of interest if applicable
if (celltype.comparison != "all"){
  sce <- sce[, sce[[celltype.col]] %in% c(celltype.of.interest, celltype.comparison)]
}

#sce$group <- sce$status

# get CPM using edgeR
SingleCellExperiment::cpm(sce) <- edgeR::cpm(counts(sce))

# get the proportion of cells with cpm > 10 for all possible combinations of cell group and gene
propOver10byGroup <- apply(
  cpm(sce),
  1,
  function(x){
  tapply(x, sce[[celltype.col]], function(x){sum(x > 10) / length(x)})
})

propOver10byGroup <- reshape2::melt(
  propOver10byGroup, varnames = c("group", "feature"), value.name = "Proportion")

# boolean for whether the proportion is > .25
propOver10byGroup$PropOver25pct <- (propOver10byGroup$Proportion > 0.25)

# now we have a data.frame with four columns - group, feature, proportion and PropOver25pct

groupsOver25pct <- as.data.frame(with(
    propOver10byGroup,
    tapply(PropOver25pct, feature, function(x){sum(x)}))
)
colnames(groupsOver25pct) <- "Groups"
groupsOver25pct$feature <- rownames(groupsOver25pct)

detectedCPMandProp <- merge(
  subset(propOver10byGroup, Proportion > 0.25), groupsOver25pct, by = "feature"
)

# detectedCPMandProp <- merge(
#   detectedCPMandProp,
#   as.data.frame(unique(colData(sce)[,c("group","time","status")])),
#   by = "group"
# )

detectedFeatures <- names(which(tapply(
  propOver10byGroup$Proportion, propOver10byGroup$feature,
  function(x){sum(x > 0.25) > 0})))

keepEndogenous <- (rownames(sce) %in% detectedFeatures)
sce <- sce[keepEndogenous]


groups <- sce[[celltype.col]]
#groups <- sce$infected
lfc <- as.numeric(snakemake@wildcards[["lfc"]])
#lfc <- 1
pval.type <- snakemake@wildcards[["pval.type"]]
#pval.type <- "all"
block <- sce[[snakemake@wildcards[["block"]]]]
#block <- sce$plate
# # calculate markers using t-test
t.test.markers <- findMarkers(sce, groups=groups, lfc=lfc, pval.type=pval.type, block=block)
t.test.df <- t.test.markers[[celltype.of.interest]]
write.table(t.test.df, file=snakemake@output[[1]], sep="\t")
# # calculate markers using wilcoxon rank sum test
wilcox.markers <- findMarkers(sce, groups=groups, test="wilcox", lfc=lfc, pval.type=pval.type, block=block)
wilcox.df <- wilcox.markers[[celltype.of.interest]]
write.table(wilcox.df, file=snakemake@output[[2]], sep="\t")
