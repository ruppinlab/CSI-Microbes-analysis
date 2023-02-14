import pandas as pd

df_3p = pd.read_csv("output/filtered-3p-units.tsv", sep="\t")
df_3p = df_3p.rename(columns={"orig.ident": "sample", "celltype": "celltype1"})
df_3p["barcode"] = df_3p.index.map(lambda x: x.split("_")[-1])
df_3p["chemistry"] = "3' v3"
df_3p = df_3p[["sample", "barcode", "celltype1", "exposure", "chemistry"]]

df_5p = pd.read_csv("output/filtered-5p-units.tsv", sep="\t")
df_5p["sample"] = "SCAF2965_5_Live"
df_5p["barcode"] = df_5p.index
df_5p["chemistry"] = "5' v2"
df_5p["exposure"] = "Live"
df_5p = df_5p.rename(columns={"celltype": "celltype1"})
df_5p = df_5p[["sample", "barcode", "celltype1", "exposure", "chemistry"]]

# desired output columns
# patient	sample	barcode	celltype1 exposure chemistry
df = pd.concat([df_3p, df_5p]).reset_index(drop=True)
df["patient"] = "P1"

df[['patient', 'sample', 'barcode', 'celltype1', 'exposure', 'chemistry']].to_csv("output/filtered-units.tsv", sep="\t", index=False)
