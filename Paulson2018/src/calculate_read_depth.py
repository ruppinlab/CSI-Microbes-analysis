import pandas as pd

meta_df = pd.read_csv("data/units.tsv", sep="\t")

# -1 means PreRxTumor and -2 means Day615Tumor
read_df = pd.read_csv("data/2586-4-tumor_exp_matrix.csv", index_col=0)
read_df = read_df.transpose()
read_df.index = read_df.index.map(lambda x: "PreRxTumor-{}".format(x) if x.endswith("-1") else "Day615Tumor-{}".format(x))
read_df.index = read_df.index.str.replace("-2", "-1")
read_df.sum().sum() # 30536526
prerx_read_df = read_df.loc[read_df.index.str.startswith("PreRxTumor")]
prerx_read_df.sum().sum() # 9161969
postrx_read_df = read_df.loc[read_df.index.str.startswith("Day615Tumor")]
postrx_read_df.sum().sum() # 21374557
# what about only tumor cells?
tumor_prerx_2586_meta_df = meta_df.loc[(meta_df["sample"] == "PreRxTumor") & (meta_df.Tumor == "Tumor")]
tumor_postrx_2586_meta_df = meta_df.loc[(meta_df["sample"] == "Day615Tumor") & (meta_df.Tumor == "Tumor")]
prerx_read_df.loc[tumor_prerx_2586_meta_df.cell].sum().sum()
postrx_read_df.loc[tumor_postrx_2586_meta_df.cell].sum().sum()


read_9245_df = pd.read_csv("data/9245-3-tumor_exp_matrix.csv", index_col=0)
read_9245_df = read_9245_df.transpose()
read_9245_df.index = read_9245_df.index.map(lambda x: "RelapseD559PBMC-{}".format(x) if x.endswith("-1") else "RelapseD565Tumor-{}".format(x))
# only focus on cells from TME
read_9245_df = read_9245_df.loc[read_9245_df.index.str.endswith("-2")]
read_9245_df.index = read_9245_df.index.str.replace("-2", "-1")
# all cells from TME read depth
read_9245_df.sum().sum() # 82855126
tumor_9245_meta_df = meta_df.loc[(meta_df["sample"] == "RelapseD565Tumor") & (meta_df.Tumor == "Tumor")]
# only tumor cells
read_9245_df.loc[tumor_9245_meta_df.cell].sum().sum() # 77548344

# comparison is 82,855,126 UMIs and 11,141 merkel polyomavirus UMIs for 5' vs
# 9,161,969 UMIs and 394 Merkel polyomavirus UMIs for pre and 21,374,557 UMIs and 691 for post 3' v2

fivep.microbial.reads <- 11141
fivep.reads <- 82855126
threep.microbial.reads.pre <- 394
threep.reads.pre <- 9161969
threep.microbial.reads.post <- 691
threep.reads.post <- 21374557

mat <- matrix(c(fivep.microbial.reads, fivep.reads, threep.microbial.reads.pre, threep.reads.pre), nrow=2, ncol=2, byrow=TRUE)
mat <- matrix(c(fivep.microbial.reads, fivep.reads, threep.microbial.reads.post, threep.reads.post), nrow=2, ncol=2, byrow=TRUE)
