import pandas as pd
#import matplotlib.pyplot as plt

patients = pd.read_csv(snakemake.input[0], sep="\t")
samples = pd.read_csv(snakemake.input[1], sep="\t")
samples = samples.merge(patients, on="patient").drop_duplicates()

#microbe_of_interest = snakemake.wildcards["microbe_of_interest"]  # "Salmonella_enterica_subsp._enterica"
output = []
for _, sample in samples.iterrows():
    try:
        pathseq_df = pd.read_csv("data/PathSeq/{}-{}/pathseq.txt".format(sample["patient"], sample["sample"]), sep="\t")
        pathseq_df = pathseq_df.loc[pathseq_df.name == "Salmonella_enterica"]

        pathseq_df["sample"] = sample["sample"]
        pathseq_df["patient"] = sample["patient"]
        output.append(pathseq_df)
    except:
        pass

df = pd.concat(output)
print(df)
df = df.loc[df.type == "species"]
df = df.sort_values(by="unambiguous", ascending=False)
df.to_csv(snakemake.output[0], sep="\t")
