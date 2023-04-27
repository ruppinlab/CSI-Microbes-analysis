library(dplyr)
library(scater)
library(scran)

calculate.Fishers.exact <- function(tax.id) {
  celltype.of.interest.positive <- dim(sce[,colData(sce)[[celltype.col]] == celltype.of.interest & counts(sce)[tax.id,] >= min.umis])[[2]]
  celltype.comparison.positive <- dim(sce[,colData(sce)[[celltype.col]] == celltype.comparison & counts(sce)[tax.id,] >= min.umis])[[2]]
  celltype.of.interest.negative <- dim(sce[,colData(sce)[[celltype.col]] == celltype.of.interest & counts(sce)[tax.id,] < min.umis])[[2]]
  celltype.comparison.negative <- dim(sce[,colData(sce)[[celltype.col]] == celltype.comparison & counts(sce)[tax.id,] < min.umis])[[2]]
  mat <- matrix(c(celltype.of.interest.positive, celltype.of.interest.negative, celltype.comparison.positive, celltype.comparison.negative), nrow=2, ncol=2, byrow=TRUE)
  f.test <- fisher.test(mat)
  return(c("p.value"=f.test$p.value[[1]], "summary.odds.ratio"=f.test$estimate[[1]], "or.ci.low"=f.test$conf.int[[1]], "or.ci.high"=f.test$conf.int[[2]], "n.reads"=sum(counts(sce)[tax.id,]), "n.positive.cells.of.interest"=celltype.of.interest.positive, "n.negative.cells.of.interest"=celltype.of.interest.negative, "n.positive.cells.comparison"=celltype.comparison.positive,
  "n.negative.cells.comparison"=celltype.comparison.negative))
}

min.umis <- as.integer(snakemake@wildcards[["min_umis"]])
counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
# counts <- read.table("output/Pt0/GSM3454529/genus_PathSeq_Bacteria_reads.tsv", sep="\t", header=TRUE, row.names=1)
# output/61_species_PathSeq_metadata.tsv
pdata <- read.table(snakemake@input[[2]], sep="\t", header = TRUE, row.names=1)
# pdata <- read.table("output/Pt0/GSM3454529/genus_PathSeq_Bacteria_metadata.tsv", sep="\t", header = TRUE, row.names=1)
row.names(pdata) <- gsub("-", ".", row.names(pdata))

column.names <- intersect(unique(colnames(counts)), unique(row.names(pdata)))
# remove samples with no spike-ins

counts <- counts[, column.names]
pdata <- pdata[column.names,]

tax.map <- read.table(snakemake@input[[3]], sep="\t", header=TRUE)
# tax.map <- read.table("output/Pt0/tax_id_map_Bacteria_PathSeq.tsv", sep="\t", header=TRUE)

celltype.col <- snakemake@wildcards[["celltype"]]
# celltype.col <- "Monocyte"
celltype.of.interest <- snakemake@wildcards[["celltype_of_interest"]]
# celltype.of.interest <- "Monocyte"
celltype.comparison <- snakemake@wildcards[["celltype_comparison"]]



# first let's remove any rows without at least five UMIs with the taxa present
counts <- counts[rowSums(counts) >= 5,]
# next, let's remove any rows without at least one entry with >= min_umis
counts <- counts[apply(counts, 1, function(x) any(x >= min.umis)),]
# if there are no reads, then write empty data.frame
if (nrow(counts) == 0){
  print("There are no reads at all. Writing empty data.frame")
  df <- data.frame(p.value=double(), summary.odds.ratio=double(), or.ci.low=double(), or.ci.high=double(), taxa=character())
  write.table(df, file=snakemake@output[[1]], sep="\t")
  quit()
}

sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData=pdata)
#print(sce)
min.proportion <- as.numeric(snakemake@wildcards[["minprop"]])

# filter down to the cell-types of interest
if (celltype.comparison != "all"){
  sce <- sce[, sce[[celltype.col]] %in% c(celltype.of.interest, celltype.comparison)]
}

# get the proportion of cells with count > 0 for all possible combinations of cell group and gene
propOver0byGroup <- apply(
  counts(sce),
  1,
  function(x){
  tapply(x, sce[[celltype.col]], function(x){sum(x > 0) / length(x)})
})
#print(propOver0byGroup)
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

#print(sce)
#print(assays(sce)$counts)
# check again after filtering
if (nrow(assays(sce)$counts) == 0){
  print("There are no reads at all. Writing empty data.frame")
  # "p.value"	"odds.ratio"	"or.ci.low"	"or.ci.high"	"taxa"
  df <- data.frame(p.value=double(), summary.odds.ratio=double(), or.ci.low=double(), or.ci.high=double(), taxa=character())
  write.table(df, file=snakemake@output[[1]], sep="\t")
  quit()
}

out <- sapply(rownames(sce), calculate.Fishers.exact)
df <- as.data.frame(t(out))
#print(df)
#print(rownames(df))
#print(tax.map)
#tax.map$tax_id <- sapply(tax.map$tax_id, toString)
df$taxa <- sapply(rownames(df), function(x) tax.map[tax.map$tax_id == x, "name"])
print(df)
write.table(df, file=snakemake@output[[1]], sep="\t")
