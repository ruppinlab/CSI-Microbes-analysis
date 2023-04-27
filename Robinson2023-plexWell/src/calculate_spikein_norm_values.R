library(scater)
library(scran)
#library(EnhancedVolcano)


source(snakemake@params[["spike_functions"]])

# sce <- generate_sce(microbe.file="output/TH179/class_PathSeq_Bacteria_reads.tsv",
# spikein.file="output/star/TH179/spike_ins.tsv", pdata.file="output/TH179/class_PathSeq_Bacteria_metadata.tsv", celltype.col="celltype1")
sce <- generate_sce(microbe.file=snakemake@input[[1]], spikein.file=snakemake@input[[2]],
  pdata.file=snakemake@input[[3]], celltype.col="Cell_Type")

tax.map <- read.table(snakemake@input[[4]], sep="\t", header=TRUE)
# tax.map <- read.table("output/TH179/tax_id_map_Bacteria_PathSeq.tsv", sep="\t", header=TRUE)
#species.tax.map <- tax.map[tax.map$taxa_level == snakemake@wildcards[["tax_level"]],]
print(sce)
# rename rownames from tax_id to names
#
# celltype.col <- snakemake@wildcards[["celltype"]]
# celltype.of.interest <- snakemake@wildcards[["celltype_of_interest"]]
# celltype.comparison <- snakemake@wildcards[["celltype_comparison"]]
#
# # now subset the cells by celltypes of interest if applicable
# if (celltype.comparison != "all"){
#   sce <- sce[, sce[[celltype.col]] %in% c(celltype.of.interest, celltype.comparison)]
# }
#
# min.proportion <- as.numeric(snakemake@wildcards[["minprop"]])
# min.cpm <- 10
#
# # get CPM using edgeR
# SingleCellExperiment::cpm(sce) <- edgeR::cpm(counts(sce))

# # if there is only one unique value in sce[[celltype.col]], this function returns a list
# # if there is > 1 unique value in sce[[celltype.col]], this function returns a matrix
# # get the proportion of cells with cpm > for all possible combinations of cell group and gene
# propOverXbyGroup <- apply(
#   cpm(sce),
#   1,
#   function(x){
#   tapply(x, sce[[celltype.col]], function(x){sum(x > min.cpm) / length(x)})
# })
#
# if (length(unique(sce[[celltype.col]])) == 1){
#   propOverXbyGroup <- t(as.matrix(propOverXbyGroup))
#   rownames(propOverXbyGroup) <- unique(sce[[celltype.col]])
# }
#
# propOverXbyGroup <- reshape2::melt(
#   propOverXbyGroup, varnames = c("group", "feature"), value.name = "Proportion")
#
# # boolean for whether the proportion is > .25
# propOverXbyGroup$PropOverXpct <- (propOverXbyGroup$Proportion > min.proportion)
#
# # now we have a data.frame with four columns - group, feature, proportion and PropOver25pct
#
# groupsOverXpct <- as.data.frame(with(
#     propOverXbyGroup,
#     tapply(PropOverXpct, feature, function(x){sum(x)}))
# )
#
# colnames(groupsOverXpct) <- "Groups"
# groupsOverXpct$feature <- rownames(groupsOverXpct)
#
# detectedFeatures <- names(which(tapply(
#   propOverXbyGroup$Proportion, propOverXbyGroup$feature,
#   function(x){sum(x > min.proportion) > 0})))
#
# keepEndogenous <- (rownames(sce) %in% detectedFeatures)
# sce <- sce[keepEndogenous]


# celltype <- sce[[celltype.col]]
# print(sce)

write.table(logcounts(sce), snakemake@output[[1]], sep="\t")
