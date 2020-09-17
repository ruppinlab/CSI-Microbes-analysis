import pandas as pd


cells = pd.read_csv(snakemake.input[0], sep="\t")
#samples = pd.read_csv(snakemake.input[1], sep="\t")
#samples = samples.merge(patients, on="patient").drop_duplicates()
# cells.set_index(["sample", "barcode"], inplace=True)
tax_level = snakemake.wildcards["tax_level"]
kingdom = snakemake.wildcards["kingdom"]
output = []

for _, cell in cells.iterrows():
    try:
        #print(cell)
        filename = snakemake.params[0].format(cell["patient"], cell["sample"], cell["barcode"])
        pathseq_df = pd.read_csv(filename, sep="\t")
        cell_name = "{}-{}".format(cell["sample"], cell["barcode"])
        pathseq_df["cell"] = cell_name
        pathseq_df = pathseq_df.loc[pathseq_df["type"] == tax_level]
        pathseq_df = pathseq_df.loc[pathseq_df["kingdom"] == kingdom]
        pathseq_df = pathseq_df.drop(columns=["tax_id", "taxonomy", "type", "kingdom", "score", "score_normalized", "reads", "reference_length"])
        if pathseq_df.empty:
            print("EMPTY")
            pathseq_df = pd.DataFrame(data={"cell": [cell_name], "name": ["placeholder"], "unambiguous": [0]})
        output.append(pathseq_df)
    except OSError as err:
        print("OS error: {0}".format(err))
        cell_name = "{}-{}".format(cell["sample"], cell["barcode"])
        print("missing {}".format(cell_name))
        cells.drop(cell.name, inplace=True)
    except:
        cell_name = "{}-{}".format(cell["sample"], cell["barcode"])
        print("missing {}".format(cell_name))
        print("Unexpected error:", sys.exc_info()[0])
        raise

df = pd.concat(output)
df = df.astype({'unambiguous': 'int64'})
# desired output - microbes are the indices and samples are the columns
read_df = df.pivot_table(index="name", columns="cell").fillna(0)
read_df.index.name = None  # remove index name
read_df.columns = read_df.columns.droplevel()  # remove multi-index
read_df.index = read_df.index.map(lambda x: x.replace("\'", "").replace("#", ""))
read_df = read_df.sort_index(axis=1)

read_df = read_df[(read_df.T != 0).any()]  # remove any rows with all zeroes
read_df.to_csv(snakemake.output[0], sep="\t")
# cells = cells.reset_index()
cells["batch"] = cells["sample"]
# use sample to refer to a single cell (which is what we do for the Smart-seq2 datasets)
cells["sample"] = cells["sample"] + "-" + cells["barcode"]
cells = cells.set_index(keys="sample")
cells = cells.sort_index()
cells.to_csv(snakemake.output[1], sep="\t")
#star_readcount_df[read_df.columns].sum().to_csv(snakemake.output[2], sep="\t")
