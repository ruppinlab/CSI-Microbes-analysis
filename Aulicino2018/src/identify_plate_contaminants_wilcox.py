import pandas as pd

output = []

for i in snakemake.input:
    df = pd.read_csv(i, sep="\t")
    df = df.loc[df["FDR"] < .05]
    output.append(df)

df = pd.concat(output)
df = df.index.str.replace("Bacteria-", "").to_frame("name")
# print(df)
df.drop_duplicates().to_csv(snakemake.output[0], sep="\t", index=False, header=False)
