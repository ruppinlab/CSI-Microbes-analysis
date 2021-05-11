import pandas as pd

tax_levels = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]

assert(len(tax_levels) == len(snakemake.input))

output = []

for i, tax_level in zip(snakemake.input, tax_levels):
    df = pd.read_csv(i, sep="\t", index_col=0)
    df["tax_level"] = tax_level
    df = df.loc[df.FDR < .05]
    output.append(df)

pd.concat(output).to_csv(snakemake.output[0], sep="\t")
