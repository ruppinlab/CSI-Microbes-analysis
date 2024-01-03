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

# get the sequencing depth for 10x 5' live exposed cells
data <- Read10X_h5("raw/P1-SCAF2965_5_Live/outs/raw_feature_bc_matrix.h5")
live <- CreateSeuratObject(counts = data, project="live", min.cells = 3, min.features=200)
units <- read.table("data/units.tsv", sep="\t", header=TRUE)
units <- units[units$sample == "SCAF2965_5_Live",]
subset_seurat_object <- subset(live, cells = units$barcode)
total.UMIs.5p <- sum(Matrix::colSums(subset_seurat_object))
# total_umis is 298,922,300 for 10x 5'
fuso.UMIs.5p <- 926
# get the sequencing depth for 10x 3' live exposed cells
data <- Read10X_h5("raw/P1-SCAF2963_3_Live/outs/raw_feature_bc_matrix.h5")
live <- CreateSeuratObject(counts = data, project="live", min.cells = 3, min.features=200)
units <- read.table("data/units.tsv", sep="\t", header=TRUE)
units <- units[units$sample == "SCAF2963_3_Live",]
subset_seurat_object <- subset(live, cells = units$barcode)
total.UMIs.3p <- sum(Matrix::colSums(subset_seurat_object))
fuso.UMIs.3p <- 1401
# total_umis is 173,100,479 for 10x 3' 

data <- matrix(c(total.UMIs.5p, fuso.UMIs.5p, total.UMIs.3p, fuso.UMIs.3p), nrow=2, byrow=TRUE)
colnames(data) <- c("human.UMIs", "Fuso.UMIs")
rownames(data) <- c("5p", "3p")

result <- chisq.test(data)