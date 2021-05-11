import pandas as pd

samples = ["BT{}".format(x) for x in range(1290, 1302)]

# cluster_df contains mapping from cell barcode to cell type annotation
cluster_df = pd.read_excel("data/6149-MetaData.xlsx", index_col=0)
# metadata_df contains mapping from sample used to download and
# sample used in cell_cluster_df
metadata_df = pd.read_csv("data/TableS2.csv")
cluster_df = cluster_df.merge(metadata_df, left_on="ID form Table S2", right_on="SAMPLE_ID")
cluster_df["barcode"] = cluster_df["cell"].apply(lambda x: x.split("_")[0])

for sample in samples:
    print(sample)
    sample_df = cluster_df.loc[cluster_df["SAMPLE_NAME"] == sample]

    try:
        df = pd.read_csv("data/{}/PathSeq_STAR_reads.tsv".format(sample), sep="\t")
    except:
        print("missing sample: {}".format(sample))
        continue

    new_df = sample_df.merge(df, left_on="barcode", right_on="CB")
    print(new_df.drop_duplicates(subset="UB"))
    print(new_df.loc[~new_df.YP.str.contains(",")].drop_duplicates(subset="UB"))


# for each different sample from a patient, we need to read in the corresponding
# tsv file containing the

# we currently identify samples using BTXXX - how do we map from BTXXXX to
# MetaData identifies samples using ID from Table S2, Patient_piece
