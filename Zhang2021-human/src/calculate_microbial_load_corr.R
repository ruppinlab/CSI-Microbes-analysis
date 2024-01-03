library(Seurat)


tumor.obj <- readRDS("output/seurat_objects/celltype1_Myeloid_cells.rds")

# subset where cells have non-zero values for the microbe of interest
subset.obj <- subset(tumor.obj, Pseudomonas > 1)


ge.mat <- as.matrix(GetAssayData(subset.obj))

cor_test_results <- apply(ge.mat, 1, function(col) cor.test(col,as.matrix(subset.obj[[genera]])[,1] , method = "spearman"))

correlations <- sapply(cor_test_results, function(result) result$estimate)
p_values <- sapply(cor_test_results, function(result) result$p.value)


df <- data.frame(pval=p_values, rho=correlations)