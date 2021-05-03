import pandas as pd

meta_df = pd.read_csv(snakemake.input[0], sep="\t")

df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
print(meta_df)
print(df)
print(snakemake.wildcards["genome"])
df = df[meta_df.loc[meta_df.infection == snakemake.wildcards["genome"], "cell"].values].sum(axis=1).sort_values(ascending=False)

df = df.loc[df >= 5]

df.to_csv(snakemake.output[0], sep="\t")
