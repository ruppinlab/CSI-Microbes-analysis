library(scater)
library(scran)


# human.counts <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
human.counts <- read.table("output/star/Pt0/human_genes.tsv", sep="\t", header=TRUE, row.names=1)
#spikein.counts <- read.table(snakemake@input[[2]], sep="\t", header=TRUE, row.names=1)
spikein.counts <- read.table("output/star/Pt0/spike_ins.tsv", sep="\t", header=TRUE, row.names=1)
microbe.counts <- read.table("output/Pt0/genus_PathSeq_Bacteria_reads.tsv", sep="\t", header=TRUE, row.names=1)
column.names <- intersect(unique(colnames(human.counts)), unique(colnames(spikein.counts)))
column.names <- intersect(column.names, unique(colnames(microbe.counts)))
# remove samples with no spike-ins
column.names <- column.names[!(column.names %in% c("B1_B003648", "P5_B000420"))]
spikein.counts <- spikein.counts[, column.names]
human.counts <- human.counts[, column.names]
microbe.counts <- microbe.counts[, column.names]
#pdata <- read.table(snakemake@input[[3]], sep="\t", header = TRUE, row.names=1)
pdata <- read.table("output/Pt0/genus_PathSeq_Bacteria_metadata.tsv", sep="\t", header = TRUE, row.names=1)
pdata <- pdata[column.names, ]
row.names(pdata) <- gsub("-", ".", row.names(pdata))

# create SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(counts = as.matrix(human.counts)), colData=pdata)
spike_se <- SummarizedExperiment(list(counts=spikein.counts))
altExp(sce, "spike") <- spike_se
microbe_se <- SummarizedExperiment(list(counts=microbe.counts))
altExp(sce, "microbe") <- microbe_se

sce <- computeSpikeFactors(sce, "spike")
sce <- logNormCounts(sce)
altExp(sce, "microbe") <- logNormCounts(altExp(sce, "microbe"), sizeFactors(sce))

# first, let's subset down to only infected and bystander cells
sce <- sce[, sce$infected %in% c("infected", "bystander")]

# now we have the spike-in normalized log read counts for all the taxa
# now, let's subset by plate
sce.p1 <- subset(sce, , plate=="P1")
median.Salmonella.p1 <- median(assays(altExp(sce.p1, "microbe"))$logcounts["590",])
sce.p1$Salmonella <- ifelse(assays(altExp(sce.p1, "microbe"))$logcounts["590",] > median.Salmonella.p1, "High", "Low")

sce.p2 <- subset(sce, , plate=="P2")
median.Salmonella.p2 <- median(assays(altExp(sce.p2, "microbe"))$logcounts["590",])
sce.p2$Salmonella <- ifelse(assays(altExp(sce.p2, "microbe"))$logcounts["590",] > median.Salmonella.p2, "High", "Low")

sce.p3 <- subset(sce, , plate=="P3")
median.Salmonella.p3 <- median(assays(altExp(sce.p3, "microbe"))$logcounts["590",])
sce.p3$Salmonella <- ifelse(assays(altExp(sce.p3, "microbe"))$logcounts["590",] > median.Salmonella.p3, "High", "Low")

sce.p4 <- subset(sce, , plate=="P4")
median.Salmonella.p4 <- median(assays(altExp(sce.p4, "microbe"))$logcounts["590",])
sce.p4$Salmonella <- ifelse(assays(altExp(sce.p4, "microbe"))$logcounts["590",] > median.Salmonella.p4, "High", "Low")

sce <- cbind(sce.p1, sce.p2, sce.p3, sce.p4)

median.Salmonella <- median(assays(altExp(sce, "microbe"))$logcounts["590",])
sce$Salmonella1 <- ifelse(assays(altExp(sce, "microbe"))$logcounts["590",] > median.Salmonella, "High", "Low")

# celltype.col <- snakemake@wildcards[["celltype"]]
celltype.col <- "Salmonella"
# celltype.of.interest <- snakemake@wildcards[["celltype_of_interest"]]
celltype.of.interest <- "High"
#celltype.comparison <- snakemake@wildcards[["celltype_comparison"]]
celltype.comparison <- "Low"

# # now subset the cells by celltypes of interest if applicable
# if (celltype.comparison != "all"){
#   sce <- sce[, sce[[celltype.col]] %in% c(celltype.of.interest, celltype.comparison)]
# }

sce$group <- sce$status

# get CPM using edgeR
SingleCellExperiment::cpm(sce) <- edgeR::cpm(counts(sce))

# get the proportion of cells with cpm > for all possible combinations of cell group and gene
propOver10byGroup <- apply(
  cpm(sce),
  1,
  function(x){
  tapply(x, sce$group, function(x){sum(x > 10) / length(x)})
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

detectedCPMandProp <- merge(
  detectedCPMandProp,
  as.data.frame(unique(colData(sce)[,c("group","time","status")])),
  by = "group"
)

detectedFeatures <- names(which(tapply(
  propOver10byGroup$Proportion, propOver10byGroup$feature,
  function(x){sum(x > 0.25) > 0})))

keepEndogenous <- (rownames(sce) %in% detectedFeatures)
sce <- sce[keepEndogenous]

# select genes with average expression level above 10 CPM in at least 25% of at least 1 group of cells
sce.2h <- sce[, sce$time == "2h"]
sce.4h <- sce[, sce$time == "4h"]
sce.6h <- sce[, sce$time == "6h"]
# print(sce)
# filterCounts <- function(m, counts = 10, cells = 10){
#   apply(m, 1, function(e){
#     return(sum(e >= counts) >= cells)
#   })
# }

#clusters <- quickCluster(sce, min.size=min(table(sce$status)), assay.type="counts")
#sce$quickCluster <- clusters
#print(sce)
#write.table(colData(sce)[, c("quickCluster", "infected")], "output/cluster_infected_4h_pi.tsv")
#print(sce[, c("quickCluster", "infected")])
#print(clusters)

# filter so we are only comparing highly expressed OTUs to reduce FDR penalty
# agg <- aggregate(colnames(sce), by=list(sce[[celltype.col]]), FUN=length)
# min.cells <- min(agg$x)*.25
# print(min.cells)
# sce <- sce[rowSums(logcounts(sce) > 5) > min.cells, ]
# # celltype.of.interest <- snakemake@wildcards[["celltype_of_interest"]]
# groups <- sce[[celltype.col]]
groups <- sce$Salmonella
# lfc <- as.numeric(snakemake@wildcards[["lfc"]])
lfc <- 1
# pval.type <- snakemake@wildcards[["pval.type"]]
pval.type <- "all"
# block <- sce[[snakemake@wildcards[["block"]]]]
block <- sce$plate
# # calculate markers using t-test
t.test.markers <- findMarkers(sce, groups=groups, lfc=lfc, pval.type=pval.type, block=block)
# t.test.df <- t.test.markers[[celltype.of.interest]]
# write.table(t.test.df, file=snakemake@output[[1]], sep="\t")
# # calculate markers using wilcoxon rank sum test
wilcox.markers <- findMarkers(sce, groups=groups, test="wilcox", lfc=lfc, pval.type=pval.type, block=block)
# wilcox.df <- wilcox.markers[[celltype.of.interest]]
# write.table(wilcox.df, file=snakemake@output[[2]], sep="\t")
