import pandas as pd


df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0) # , nrows=100
clinical_df = pd.read_excel(snakemake.input[1], sheet_name="ExtraEndpoints")

#clinical_df = clinical_df.loc[clinical_df["type"].isin(["LUAD", "LUSC"])]
clinical_df = clinical_df.loc[clinical_df["type"] == snakemake.wildcards["cancertype"]]

# select only LUAD and LUSC patients
df = df[df.columns[df.columns.map(lambda x: x.startswith(tuple(clinical_df["bcr_patient_barcode"])))]]
# select only one sample per patient (and make sure it is a tumor sample!)
df = df.T
df["sample"] = df.index.map(lambda x: x.split("-")[3][0:2])
# df.columns = df.columns.map(lambda x: "-".join(x.split("-")[0:3]))
df = df.astype({"sample": "int32"})
# technically < 10 is a tumor sample but we only want one sample per tumor
# and there are only two samples where sample == 2 and each of them have another sample == 1
df = df.loc[df["sample"] == 1]
df = df.T
df["gene"] = df.index.map(lambda x: x.split("|")[0])
df = df.loc[~df["gene"].duplicated(keep=False)]
df = df.set_index("gene")
df = df.drop(index="sample")
df.to_csv(snakemake.output[0], sep="\t")
