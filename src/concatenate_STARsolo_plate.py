import pandas as pd

output = []
for i in snakemake.input:
    df = pd.read_csv(i, index_col=0, sep="\t")
    output.append(df)

df = pd.concat(output, axis=1)

df.to_csv(snakemake.output[0], sep="\t")
