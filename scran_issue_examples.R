
# example #1 - this does not generate significant differences - works as expected

counts <- c(rep(5, 10), rep(20, 10), rep(20, 10), rep(5, 10))
mat <- as.matrix(t(counts))
rownames(mat) <- "g1"
#colnames(mat) <- seq(1, 20)
#samples <- seq(1, 20)
batch <- c(rep("B1", 20), rep("B2", 20))
group <- c(rep("G1", 10), rep("G2", 10), rep("G1", 10), rep("G2", 10))
df <- data.frame(batch=batch, group=group)
sce <- SingleCellExperiment(assays = list(counts = mat), colData=df)
out <- findMarkers(sce, group=sce$group, block=sce$block, lfc=0.5, direction="any", assay.type = "counts")

# example #2 - should generate significant DE - it does
counts <- c(rep(5, 10), rep(20, 10), rep(5, 10), rep(20, 10))
mat <- as.matrix(t(counts))
rownames(mat) <- "g1"
#colnames(mat) <- seq(1, 20)
#samples <- seq(1, 20)
batch <- c(rep("B1", 20), rep("B2", 20))
group <- c(rep("G1", 10), rep("G2", 10), rep("G1", 10), rep("G2", 10))
df <- data.frame(batch=batch, group=group)
sce <- SingleCellExperiment(assays = list(counts = mat), colData=df)
out <- findMarkers(sce, group=sce$group, block=sce$block, lfc=0.5, direction="any", assay.type = "counts")

# example 3
counts <- c(rep(5, 10), rep(20, 10), rep(12, 10), rep(10, 10))
mat <- as.matrix(t(counts))
rownames(mat) <- "g1"
#colnames(mat) <- seq(1, 20)
#samples <- seq(1, 20)
batch <- c(rep("B1", 20), rep("B2", 20))
group <- c(rep("G1", 10), rep("G2", 10), rep("G1", 10), rep("G2", 10))
df <- data.frame(batch=batch, group=group)
sce <- SingleCellExperiment(assays = list(counts = mat), colData=df)
out <- findMarkers(sce, group=sce$group, block=sce$block, lfc=0.5, direction="any", assay.type = "counts")


# example 4
counts <- c(rep(5, 10), rep(7, 10), rep(0, 10), rep(0, 10))
mat <- as.matrix(t(counts))
rownames(mat) <- "g1"
#colnames(mat) <- seq(1, 20)
#samples <- seq(1, 20)
batch <- c(rep("B1", 20), rep("B2", 20))
group <- c(rep("G1", 10), rep("G2", 10), rep("G1", 10), rep("G2", 10))
df <- data.frame(batch=batch, group=group)
sce <- SingleCellExperiment(assays = list(counts = mat), colData=df)
out <- findMarkers(sce, group=sce$group, block=sce$block, lfc=0.5, direction="any", assay.type = "counts")
