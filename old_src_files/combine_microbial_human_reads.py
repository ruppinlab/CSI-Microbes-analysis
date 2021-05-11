import pandas as pd

human_read_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
microbe_read_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
microbe_read_df.index = microbe_read_df.index.map(lambda x: "{}-{}".format(snakemake.wildcards["kingdom"], x))

pd.concat([human_read_df[microbe_read_df.columns], microbe_read_df]).to_csv(snakemake.output[0], sep="\t")
