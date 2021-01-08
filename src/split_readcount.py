import pandas as pd

star_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)

spikein_df = star_df.loc[star_df.index.str.startswith(snakemake.params["spike"])]
human_df = star_df.loc[~star_df.index.str.startswith(snakemake.params["spike"])]
spikein_df.to_csv(snakemake.output[0], sep="\t")
human_df.to_csv(snakemake.output[1], sep="\t")
