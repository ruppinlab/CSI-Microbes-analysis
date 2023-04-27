library(Seurat)


live.read.matrix <- Read10X_h5("data/SCAF2963_3_Live_filtered_feature_bc_matrix.h5")
colnames(live.read.matrix) <- lapply(colnames(live.read.matrix), function(x){paste0("SCAF2963_3_Live-", x)})

meta.data <- read.table("output/all_celltypes_DE_metadata.tsv", sep="\t", header=TRUE, row.names = 1)

read.matrix <- live.read.matrix
cell.names <- intersect(unique(row.names(meta.data)), unique(colnames(read.matrix)))
meta.data <- meta.data[cell.names,]
read.matrix <- read.matrix[, cell.names]
sum(read.matrix) # 204510528

microbial.reads.plate <- 481184
human.reads.plate <- 301507801
microbial.reads.droplet3 <- 1401
human.reads.droplet3 <- 204510528
mat <- matrix(c(microbial.reads.plate, human.reads.plate, microbial.reads.droplet3, human.reads.droplet3), nrow=2, ncol=2, byrow=TRUE)
