import pandas as pd
import matplotlib.pyplot as plt


output = []
for input in snakemake.input:
    df = pd.read_csv(input, sep="\t", index_col=0)
    output.append(df)

df = pd.concat(output)

plt.plot(df.index, df["sig_percent"])
plt.ylabel("Percentage of samples with {} found to be DA".format(snakemake.wildcards["microbe"]))
plt.xlabel("Number of reads from {}".format(snakemake.wildcards["microbe"]))
plt.savefig(snakemake.output[0])
