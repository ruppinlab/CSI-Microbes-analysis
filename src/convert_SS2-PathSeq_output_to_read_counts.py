import pandas as pd


cells = pd.read_csv(snakemake.input[0], sep="\t", dtype={"patient": "str", "sample": "str"})
cells.set_index("cell", inplace=True)
cells = cells.loc[cells.patient == snakemake.wildcards["patient"]]
tax_level = snakemake.wildcards["tax_level"]
kingdom = snakemake.wildcards["kingdom"]
output = []
print(cells.shape)
for _, cell in cells.iterrows():
    try:

        filename = snakemake.params[0].format(cell.patient, cell["sample"], cell.plate, cell.name)
        pathseq_df = pd.read_csv(filename, sep="\t")
        pathseq_df["sample"] = cell.name
        pathseq_df = pathseq_df.loc[pathseq_df["type"] == tax_level]
        if kingdom != "All":
            if kingdom == "Viruses":
                pathseq_df = pathseq_df.loc[pathseq_df["taxonomy"].str.startswith("root|Viruses")]
            else:
                pathseq_df = pathseq_df.loc[pathseq_df["kingdom"] == kingdom]
        pathseq_df = pathseq_df.drop(columns=["name", "taxonomy", "type", "kingdom", "score", "score_normalized", "reads", "reference_length"])
        pathseq_df = pathseq_df.rename(columns={"tax_id": "name"})
        if pathseq_df.empty:
            pathseq_df = pd.DataFrame(data={"sample": [cell.name], "name": ["placeholder"], "unambiguous": [0]})
        output.append(pathseq_df)
    except OSError as err:
        print("OS error: {0}".format(err))
        print("missing {}".format(cell.name))
        cells.drop(cell.name, inplace=True)
    except:
        print("missing {}".format(cell.name))
        print("Unexpected error:", sys.exc_info()[0])
        raise

df = pd.concat(output)
df = df.astype({'unambiguous': 'int64'})
# desired output - microbes are the indices and samples are the columns
read_df = df.pivot_table(index="name", columns="sample").fillna(0)
read_df.index.name = None  # remove index name
read_df.columns = read_df.columns.droplevel()  # remove multi-index
#read_df.index = read_df.index.map(lambda x: x.replace("\'", "").replace("#", ""))
read_df = read_df.sort_index(axis=1)

read_df = read_df[(read_df.T != 0).any()]  # remove any rows with all zeroes
read_df.to_csv(snakemake.output[0], sep="\t")
cells = cells.sort_index()
cells.to_csv(snakemake.output[1], sep="\t")
#star_readcount_df[read_df.columns].sum().to_csv(snakemake.output[2], sep="\t")
