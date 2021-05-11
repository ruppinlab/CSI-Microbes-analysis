library(scran)
library(scater)
# library(EnhancedVolcano)

print(snakemake@wildcards[["seed"]])
set.seed(strtoi(snakemake@wildcards[["seed"]]))

sce <- readRDS(snakemake@input[[1]])

celltype.col <- "patient"
min.proportion <- .5
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

groups <- ifelse(sce$patient %in% c("TH231", "TH236", "TH238", "TH266"), "infected", "control")
# we have 4 infected tumors and 9 uninfected tumors - let's randomly assign a batch to each
infected.batch <- seq(1, 4)
uninfected.batch <- rep_len(infected.batch, 9)
infected.batch <- sample(infected.batch)
names(infected.batch) <- c("TH231", "TH236", "TH238", "TH266")
uninfected.batch <- sample(uninfected.batch)
names(uninfected.batch) <- c("TH067", "TH103", "TH171", "TH179", "TH220", "TH225", "TH226", "TH227", "TH248")
batch <- c(infected.batch, uninfected.batch)
sce$batch <- apply(colData(sce), 1, function(x) batch[[x["patient"]]])

lfc <- as.numeric(snakemake@wildcards[["lfc"]])
#pval.type <- snakemake@wildcards[["pval.type"]]
#block <- sce[[snakemake@wildcards[["block"]]]]
direction <- snakemake@wildcards[["direction"]]
# calculate markers using t-test
t.test.markers <- findMarkers(sce, groups=groups, lfc=lfc, block=sce$batch, direction=direction)
t.test.df <- t.test.markers[["infected"]]
write.table(t.test.df, file=snakemake@output[[1]], sep="\t")
# calculate markers using wilcoxon rank sum test
wilcox.markers <- findMarkers(sce, groups=groups, test="wilcox", lfc=lfc, block=sce$batch, direction=direction)
wilcox.df <- wilcox.markers[["infected"]]
write.table(wilcox.df, file=snakemake@output[[2]], sep="\t")

# EnhancedVolcano(df, lab=df$taxa, x = "summary.logFC", y = "p.value",
#   FCcutoff = .5, pCutoff = .05, selectLab = microbes.of.interest, subtitle=NULL,
#   caption=NULL, title="Infected versus Bystander", xlim=c(-2,3),
#   ylim = c(0, -log10(10e-13))
# )

#EnhancedVolcano(t.test.df, lab=rownames(t.test.df), x = "summary.logFC", y = "p.value", subtitle=NULL, caption=NULL, title="Infected versus uninfected tumor cells")
#ggsave(snakemake@output[[3]])
