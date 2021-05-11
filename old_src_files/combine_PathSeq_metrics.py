import pandas as pd

samples = pd.read_csv(snakemake.input[0], sep="\t")
samples = samples.loc[samples["patient"] == "BC06"]
output_df = []
for _, samp in samples.iterrows():
    try:
        fname = snakemake.params[0].format(samp["patient"], samp["sample"])
        df = pd.read_csv(fname, sep="\t", comment="#")
        df.index = [samp["sample"]]
        output_df.append(df)
    except OSError as err:
        print("OS error: {0}".format(err))
        print("missing {}".format(samp["sample"]))

pd.concat(output_df).to_csv(snakemake.output[0], sep="\t")
