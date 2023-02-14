library(scater)
library(scuttle)
library(scran)
library(ggplot2)
library(dynamicTreeCut)
library(stringr)


set.seed(1010)
# taken mostly from http://bioconductor.org/books/3.16/OSCA.workflows/lun-416b-cell-line-smart-seq2.html#lun-416b-cell-line-smart-seq2

counts <- read.table("output/star/P1/matrix.tsv")
spike.mat <- counts[grepl("^ERCC-", rownames(counts)),]
mat <- counts[!grepl("^ERCC-", rownames(counts)),]

meta <- read.table("data/unfiltered-units.tsv", header=TRUE, sep="\t")

# sort so cells are in the same order
meta <- meta[order(meta$cell),]
mat <- mat[, order(colnames(mat))]
spike.mat <- spike.mat[, order(colnames(spike.mat))]

sce <- SingleCellExperiment(list(counts=mat), colData=DataFrame(label=meta))
altExp(sce, "ERCC") <- SummarizedExperiment(list(counts=spike.mat))

# remove "label." prefix for column names
colnames(colData(sce)) <- str_replace(colnames(colData(sce)), "label.", "")

unfiltered <- sce # make a copy before filtering

# calculate and plot filtering criteria
mito <- which(startsWith(rownames(unfiltered),"MT-"))

stats <- perCellQCMetrics(unfiltered, subsets=list(Mito=mito))
qc <- perCellQCFilters(stats, sub.fields=c("subsets_Mito_percent", "altexps_ERCC_percent"))

colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$celltype <- factor(colData(unfiltered)$Cell_Type)
unfiltered$discard <- qc$discard


g <- gridExtra::grid.arrange(
    plotColData(unfiltered, x="celltype", y="sum",
        colour_by="discard") + scale_y_log10() + ggtitle("Total count"),
    plotColData(unfiltered, x="celltype", y="detected",
        colour_by="discard") + scale_y_log10() + ggtitle("Detected features"),
    plotColData(unfiltered, x="celltype", y="subsets_Mito_percent",
        colour_by="discard") + ggtitle("Mito percent"),
    plotColData(unfiltered, x="celltype", y="altexps_ERCC_percent",
        colour_by="discard") + ggtitle("ERCC percent"),
    nrow=1,
    ncol=4
)
ggsave("output/plots/human_QC_filtering_plot.pdf", g, width=20, height=4)

# to visulize why cells were removed
colSums(as.matrix(qc))
# now actually filter
sce$discard <- qc$discard
sce <- sce[, !sce$discard]
# let's also discard CT-HCT116 since we are not using for downstream analysis
sce <- sce[,colData(sce)$Cell_Type %in% c("Jurkat", "HCT116")]
# also drop the cells with no ERCCs
sce <- sce[,colData(sce)$Cell_Type_and_Condition != "noERCC-HCT116-infected"]
# be sure to include the empty wells in the output units file
empty.wells <- unfiltered[,unfiltered$Cell_Type == "ERCC only"]
print(empty.wells)
units1 <- colData(sce)[c("patient", "sample", "plate", "cell", "Cell_Type", "Cell_Type_and_Condition", "Condition")]
units2 <- colData(empty.wells)[c("patient", "sample", "plate", "cell", "Cell_Type", "Cell_Type_and_Condition", "Condition")]
units <- rbind(units1, units2)
write.table(units, "output/filtered-units.tsv", sep="\t")
# perform normalization
sce <- computeSpikeFactors(sce, "ERCC")
sce <- logNormCounts(sce)

# variance modeling
dec <- modelGeneVarWithSpikes(sce, "ERCC") # no block - block=sce.416b$block
chosen.hvgs <- getTopHVGs(dec, prop=0.1)

# dimensionality reduction
sce <- runPCA(sce, ncomponents=10, subset_row=chosen.hvgs, BSPARAM=BiocSingular::ExactParam())

sce <- runTSNE(sce, dimred="PCA", perplexity=10)
sce <- runUMAP(sce, dimred="PCA")
# clustering
my.dist <- dist(reducedDim(sce, "PCA"))
my.tree <- hclust(my.dist, method="ward.D2")

my.clusters <- unname(cutreeDynamic(my.tree, distM=as.matrix(my.dist),
    minClusterSize=10, verbose=0))
colLabels(sce) <- factor(my.clusters)

saveRDS(sce, "output/sce_filtered_obj.rds")
