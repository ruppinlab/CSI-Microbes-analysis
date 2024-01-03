library(Seurat)

data <- Read10X_h5("raw/P1-SCAF2961_1_Uninfected/outs/raw_feature_bc_matrix.h5")
uninf.df <- data.frame(patient="P1", sample="SCAF2961_1_Uninfected", barcode=colnames(data))

data <- Read10X_h5("raw/P1-SCAF2962_2_HK/outs/raw_feature_bc_matrix.h5")
HK.df <- data.frame(patient="P1", sample="SCAF2962_2_HK", barcode=colnames(data))

data <- Read10X_h5("raw/P1-SCAF2963_3_Live/outs/raw_feature_bc_matrix.h5")
live.3p.df <- data.frame(patient="P1", sample="SCAF2963_3_Live", barcode=colnames(data))

data <- Read10X_h5("raw/P1-SCAF2965_5_Live/outs/raw_feature_bc_matrix.h5")
live.5p.df <- data.frame(patient="P1", sample="SCAF2965_5_Live", barcode=colnames(data))

df <- rbind(uninf.df, HK.df, live.3p.df, live.5p.df)
write.csv(df, "output/unfiltered-units.tsv", sep="\t")