import pandas as pd

output_df = []
for f in snakemake.input:
    df = pd.read_csv(f, sep="\t", index_col=0)["FDR"]
    output_df.append(df)

print(output_df)
df = pd.concat(output_df, axis=1)
print(df)
df.max(axis=1, skipna=False).sort_values().to_frame(name="FDR").to_csv(snakemake.output[0], sep="\t")
