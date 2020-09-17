import pandas as pd


patients = pd.read_csv(snakemake.input[0], sep="\t")
samples = pd.read_csv(snakemake.input[1], sep="\t")
samples = samples.merge(patients, on="patient").drop_duplicates()
samples.set_index("sample", inplace=True)
tax_level = snakemake.wildcards["tax_level"]
kingdom = snakemake.wildcards["kingdom"]
output = []
# samples = samples.iloc[0:100]
for _, sample in samples.iterrows():
    try:
        filename = snakemake.params[0].format(sample.plate, sample.name)
        pathseq_df = pd.read_csv(filename, sep="\t")
        pathseq_df["sample"] = sample.name
        pathseq_df = pathseq_df.loc[pathseq_df["type"] == tax_level]
        pathseq_df = pathseq_df.loc[pathseq_df["kingdom"] == kingdom]
        pathseq_df = pathseq_df.drop(columns=["tax_id", "taxonomy", "type", "kingdom", "score", "score_normalized", "reads", "reference_length"])
        if pathseq_df.empty:
            pathseq_df = pd.DataFrame(data={"sample": [sample.name], "name": ["placeholder"], "unambiguous": [0]})
        output.append(pathseq_df)
    except OSError as err:
        print("OS error: {0}".format(err))
        print("missing {}".format(sample.name))
        samples.drop(sample.name, inplace=True)
    except:
        print("missing {}".format(sample.name))
        print("Unexpected error:", sys.exc_info()[0])
        raise

df = pd.concat(output)
df = df.astype({'unambiguous': 'int64'})
# desired output - microbes are the indices and samples are the columns
read_df = df.pivot_table(index="name", columns="sample").fillna(0)
#print(read_df)
read_df.index.name = None  # remove index name
read_df.columns = read_df.columns.droplevel()  # remove multi-index
read_df.index = read_df.index.map(lambda x: x.replace("\'", "").replace("#", ""))
read_df = read_df.sort_index(axis=1)

read_df = read_df[(read_df.T != 0).any()]  # remove any rows with all zeroes
read_df.to_csv(snakemake.output[0], sep="\t")
samples = samples.sort_index()
samples.to_csv(snakemake.output[1], sep="\t")
#star_readcount_df[read_df.columns].sum().to_csv(snakemake.output[2], sep="\t")
