import pandas as pd


cells = pd.read_csv(snakemake.input[0], sep="\t", dtype={"patient": "str"})
cells = cells.loc[cells.patient == snakemake.wildcards["patient"]]

kingdom = snakemake.wildcards["kingdom"]
output = []
print(cells)
cells["cell"] = cells.apply(lambda x: "{}-{}".format(x["sample"], x["barcode"]), axis=1)

for _, cell in cells.iterrows():
    try:
        filename = snakemake.params[0].format(cell["patient"], cell["sample"], cell["barcode"])
        pathseq_df = pd.read_csv(filename, sep="\t")
        pathseq_df["sample"] = cell.name
        pathseq_df = pathseq_df.loc[pathseq_df["type"].isin(["superkingdom", "phylum", "class", "order", "family", "genus", "species"])]
        if kingdom != "All":
            if kingdom == "Viruses":
                pathseq_df = pathseq_df.loc[pathseq_df["taxonomy"].str.startswith("root|Viruses")]
            else:
                pathseq_df = pathseq_df.loc[pathseq_df["kingdom"] == kingdom]
        # get tax_id instead of name
        output.append(pathseq_df[["name", "tax_id", "type"]])
    except OSError as err:
        print("OS error: {0}".format(err))
        print("missing {}".format(cell.name))
        cells.drop(cell.name, inplace=True)
    except:
        print("missing {}".format(cell.name))
        print("Unexpected error:", sys.exc_info()[0])
        raise

df = pd.concat(output)
df["name"] = df["name"].map(lambda x: x.replace("\'", "").replace("#", ""))
#print(df)
df = df.drop_duplicates()
df = df.rename(columns={"type": "taxa_level"})
df.to_csv(snakemake.output[0], sep="\t", index=False)
#star_readcount_df[read_df.columns].sum().to_csv(snakemake.output[2], sep="\t")
