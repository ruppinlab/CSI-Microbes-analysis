import pandas as pd
import numpy as np

df = pd.read_csv(snakemake.input[0], sep=" ", index_col=0)

df["metric"] = df.apply(lambda x: np.log10(x["p.value"]) * (-1 if x["rho"] > 0 else 1), axis=1)
df = df.set_index("gene2")
df["metric"].sort_values().to_csv(snakemake.output[0], sep="\t", header=False)
