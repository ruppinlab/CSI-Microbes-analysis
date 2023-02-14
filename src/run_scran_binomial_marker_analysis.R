library(dplyr)
library(scater)
library(scran)
#library(ggplot2)
#library(EnhancedVolcano)

# counts <- read.table("output/P1/SCAF2965_5_Live/genus_PathSeq_All_reads.tsv", sep="\t", header=TRUE, row.names=1)
counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
# pdata <- read.table("output/P1/PathSeq_metadata.tsv", sep="\t", header = TRUE, row.names=1)
pdata <- read.table(snakemake@input[[2]], sep="\t", header = TRUE, row.names=1)
row.names(pdata) <- gsub("-", ".", row.names(pdata))
# tax.map <- read.table("output/P1/tax_id_map_All_PathSeq.tsv", sep="\t", header=TRUE)
tax.map <- read.table(snakemake@input[[3]], sep="\t", header=TRUE)
# min.umis <- 2
min.umis <- as.numeric(snakemake@wildcards[["min_umis"]])
column.names <- intersect(unique(colnames(counts)), unique(row.names(pdata)))
# remove samples with no spike-ins

counts <- counts[, column.names]
pdata <- pdata[column.names,]
rownames(counts) <- sapply(rownames(counts), function(x) tax.map[tax.map$tax_id == x, "name"])

# first let's remove any rows without at least five UMIs with the taxa present
counts <- counts[rowSums(counts) >= 5,]
# next, let's remove any rows without at least one entry with >= min_umis
counts <- counts[apply(counts, 1, function(x) sum(x >= min.umis) >= 5),]

celltype.col <- snakemake@wildcards[["celltype"]]
#celltype.of.interest <- snakemake@wildcards[["celltype_of_interest"]]
#celltype.comparison <- snakemake@wildcards[["celltype_comparison"]]

# if there are no reads, then write empty data.frame
if (nrow(counts) == 0){
  print("There are no reads at all. Writing empty data.frame")
  df <- data.frame(Top=character(), p.value=double(), FDR=double(), summary.logFC=double())
  # if (celltype.comparison != "all"){
  #   n <- paste0("LogFC.", celltype.comparison)
  #   df[n] = double()
  # } else {
  # cols <- setdiff(unique(pdata[[celltype.col]]), celltype.of.interest)
  # for (c in cols)
  # {
  #   n <- paste0("LogFC.", c)
  #   df[n] <- double()
  # }
  # }
  write.table(df, file=snakemake@output[[1]], sep="\t")
  quit()
}

sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData=pdata)

# now subset the cells by celltypes of interest if applicable
# if (celltype.comparison != "all"){
#   sce <- sce[, sce[[celltype.col]] %in% c(celltype.of.interest, celltype.comparison)]
# }


celltype <- sce[[celltype.col]]

#groups <- celltype
pval.type <- snakemake@wildcards[["pvaltype"]]
direction <- snakemake@wildcards[["direction"]]

binom <- findMarkers(sce, test="binom", groups=celltype, pval.type=pval.type, threshold=min.umis, assay.type="counts", direction=direction)
# binom is a list
# save the row.names as taxa
for (celltype in names(binom)){
  binom[[celltype]]["taxa"] = rownames(binom[[celltype]])
}
# convert to data.frame
binom <- lapply(binom, as.data.frame)
binom <- bind_rows(binom, .id = "celltype")

# I want to store the rownames in a column
#df <- binom[[celltype.of.interest]]
# df$taxa <- lapply(rownames(df), function(x) tax.map[tax.map$tax_id == x, "name"])
write.table(binom, file=snakemake@output[[1]], sep="\t", row.names=FALSE)
