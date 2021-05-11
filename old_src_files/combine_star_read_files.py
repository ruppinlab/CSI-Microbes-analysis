import pandas as pd


#patients = pd.read_csv(snakemake.input[0], sep="\t")
samples = pd.read_csv(snakemake.input[1], sep="\t")
#samples = samples.merge(patients, on="patient").drop_duplicates()
output = []
for _, sample in samples.iterrows():
    try:
        star_df = pd.read_csv(snakemake.params[0].format(sample["patient"], sample["sample"]), sep="\t", skiprows=4, header=None, index_col=0)
        output.append(star_df[1].rename(sample["sample"])) # unstranded
    except:
        print("missing " + snakemake.params[0].format(sample["patient"], sample["sample"]))

df = pd.concat(output, axis=1)
df.to_csv(snakemake.output[0], sep="\t")
