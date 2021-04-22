library(scran)
library(scater)
library(EnhancedVolcano)


sce <- readRDS(snakemake@input[[1]])

# filter to remove tumors with < 15 tumor cells (but keep TH266)
#sce <- sce[,sce$patient %in% c("TH067", "TH171", "TH179", "TH220", "TH226", "TH231", "TH236", "TH238", "TH248", "TH266")]

# we want to compare one infected tumor vs. all of the uninfected tumors
uninfected.tumors <- c("TH067", "TH171", "TH179", "TH220", "TH226", "TH248")
tumor.of.interest <- snakemake@wildcards[["tumor_of_interest"]]
#tumor.of.interest <- "TH236"
#tumors.of.interest <- c(uninfected.tumors, tumor.of.interest)
sce <- sce[, sce$patient %in% c(uninfected.tumors, tumor.of.interest)]

celltype.col <- "patient"
min.proportion <- .50
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

#groups <- ifelse(sce$patient %in% c("TH231", "TH236", "TH238", "TH266"), "infected", "control")
#lfc <- 0.5
lfc <- as.numeric(snakemake@wildcards[["lfc"]])
#pval.type <- "all"
pval.type <- snakemake@wildcards[["pvaltype"]]
print(pval.type)
direction <- snakemake@wildcards[["direction"]]
#block <- sce[[snakemake@wildcards[["block"]]]]
#direction <- snakemake@wildcards[["direction"]]
# calculate markers using t-test
t.test.markers <- findMarkers(sce, groups=sce$patient, lfc=lfc, pval.type=pval.type, direction=direction)
t.test.df <- t.test.markers[[tumor.of.interest]]
write.table(t.test.df, file=snakemake@output[[1]], sep="\t")
# calculate markers using wilcoxon rank sum test
wilcox.markers <- findMarkers(sce, groups=sce$patient, test="wilcox", lfc=lfc, pval.type=pval.type, direction=direction)
wilcox.df <- wilcox.markers[[tumor.of.interest]]
write.table(wilcox.df, file=snakemake@output[[2]], sep="\t")
#EnhancedVolcano(t.test.df, lab=rownames(t.test.df), x = "summary.logFC", y = "p.value")
#ggsave(snakemake@output[[3]])
