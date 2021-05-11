import pandas as pd

units = pd.read_csv(snakemake.input[0], sep="\t")
output = []
for _, sample in units.iterrows():
    try:
        star_df = pd.read_csv(snakemake.params[0].format(sample["patient"], sample["sample"], sample["plate"], sample["cell"]), sep="\t", skiprows=4, header=None, index_col=0)
        output.append(star_df[1].rename(sample["cell"])) # unstranded
    except:
        print("missing " + snakemake.params[0].format(sample["patient"], sample["sample"], sample["plate"], sample["cell"]))

df = pd.concat(output, axis=1)
df.to_csv(snakemake.output[0], sep="\t")
