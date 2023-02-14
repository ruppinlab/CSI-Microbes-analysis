import pandas as pd


cells = pd.read_csv(snakemake.input[0], sep="\t", dtype={"patient": "str"})
#print(cells)
cells = cells.loc[cells["patient"] == snakemake.wildcards["patient"]]
#samples = pd.read_csv(snakemake.input[1], sep="\t")
#samples = samples.merge(patients, on="patient").drop_duplicates()
# cells.set_index(["sample", "barcode"], inplace=True)

# tax_level = snakemake.wildcards["tax_level"]
# kingdom = snakemake.wildcards["kingdom"]
output = []
#print(cells)
cells["cell"] = cells.apply(lambda x: "{}-{}".format(x["sample"], x["barcode"]), axis=1)
#print(cells)
for _, cell in cells.iterrows():
    try:
        filename = snakemake.params[0].format(cell["patient"], cell["sample"], cell["barcode"])
        pathseq_df = pd.read_csv(filename, sep="\t")
        cell_name = cell["cell"]
        pathseq_df["cell"] = cell_name
        # pathseq_df = pathseq_df.loc[pathseq_df["type"] == tax_level]
        # if kingdom != "All":
        #     if kingdom == "Viruses":
        #         pathseq_df = pathseq_df.loc[pathseq_df["taxonomy"].str.startswith("root|Viruses")]
        #     else:
        #         pathseq_df = pathseq_df.loc[pathseq_df["kingdom"] == kingdom]
        pathseq_df = pathseq_df.drop(columns=["name", "taxonomy", "score", "score_normalized", "reads", "reference_length"])
        pathseq_df = pathseq_df.rename(columns={"tax_id": "name"})
        if pathseq_df.empty:
            # print("EMPTY")
            pathseq_df = pd.DataFrame(data={"cell": [cell_name], "name": ["placeholder"], "unambiguous": [0], "type": ["placeholder"], "kingdom": ["placeholder"]})
        # else:
        #     print(pathseq_df)
        # print(pathseq_df)
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
read_df = df.drop(columns=["type", "kingdom"]).pivot_table(index="name", columns="cell").fillna(0)
read_df = read_df[(read_df.T != 0).any()]
read_df.index.name = None  # remove index name
read_df.columns = read_df.columns.droplevel()  # remove multi-index

tax_levels = ["root", "superkingdom", "kingdom", "phylum", "class", "order",
              "family", "genus", "species", "strain", "no_rank"]

for tax_level in tax_levels:
    taxa = df.loc[(df.type == tax_level) & (df["unambiguous"] > 0), "name"].unique()
    tax_read_df = read_df.loc[taxa]
    tax_read_df = tax_read_df.sort_index(axis=1)
    out_file = snakemake.params[1].format(tax_level=tax_level, patient=snakemake.wildcards["patient"], kingdom="All", method="PathSeq")
    tax_read_df.to_csv(out_file, sep="\t")


cells = cells.set_index(keys="cell")
cells = cells.sort_index()
cells.to_csv(snakemake.output[0], sep="\t")
