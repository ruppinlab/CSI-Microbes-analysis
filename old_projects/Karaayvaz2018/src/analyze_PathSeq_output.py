import pandas as pd
#import matplotlib.pyplot as plt

patients = pd.read_csv(snakemake.input[0], sep="\t")
samples = pd.read_csv(snakemake.input[1], sep="\t")
samples = samples.merge(patients, on="patient").drop_duplicates()

output = []
for _, sample in samples.iterrows():
    try:
        pathseq_df = pd.read_csv("data/PathSeq/{}-{}/pathseq.txt".format(sample["patient"], sample["sample"]), sep="\t")

        #pathseq_df = pathseq_df.loc[pathseq_df.name == microbe_of_interest]
        pathseq_df["sample"] = sample["sample"]
        print(pathseq_df)
        #pathseq_df["status"] = sample["status"]
        #pathseq_df["infection"] = sample["infection"]
        output.append(pathseq_df)
    except:
        print("missing {}".format(sample["sample"]))

df = pd.concat(output)
df = df.loc[df.type == "species"]
df = df.loc[df.unambiguous > 100]
df = df.sort_values(by="unambiguous")
df.to_csv(snakemake.output[0], sep="\t")
# fig1, ax1 = plt.subplots(figsize=(20, 10))
# ax1.boxplot([df.loc[df.infection == "Mock", "unambiguous"],
#              df.loc[(df.infection == "D23580") & (df.status == "infected"), "unambiguous"],
#              df.loc[(df.infection == "D23580") & (df.status == "exposed"), "unambiguous"],
#              df.loc[(df.infection == "LT2") & (df.status == "infected"), "unambiguous"],
#              df.loc[(df.infection == "LT2") & (df.status == "exposed"), "unambiguous"]],
#             labels=["Mock", "D23580 infected", "D23580 exposed", "LT2 infected", "LT2 exposed"])
# ax1.set_xlabel("condition")
# ax1.set_ylabel("number of unambiguous reads per cell")
# ax1.set_title("Number of reads from {} across conditions".format(microbe_of_interest))
# #plt.show()
# plt.savefig(snakemake.output[0])
# print(df.groupby(["status", "infection"])["unambiguous"].mean())
