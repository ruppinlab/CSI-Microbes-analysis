import pandas as pd
import numpy as np

parent_taxa = {
    "species": "genus",
    "genus": "family",
    "family": "order",
    "order": "class",
    "class": "phylum",
    "phylum": "superkingdom"
}

def get_parent(x):
    #print(x)
    if x["type"] == "superkingdom":
        return np.nan
    else:
        #print(x)
        #print(x["taxonomy"])
        ancestor_name = x["taxonomy"].split("|")[-2]
        ancestor = pathseq_df.loc[pathseq_df["name"] == ancestor_name]
        #print(ancestor)
        # this means the ancestor is no_rank
        if ancestor.empty:
            ancestor_name = x["taxonomy"].split("|")[-3]
            #print(ancestor_name)
            ancestor = pathseq_df.loc[pathseq_df["name"] == ancestor_name]
            # this means the ancestor is species_group:
            if ancestor.empty:
                ancestor_name = x["taxonomy"].split("|")[-4]
                #print(ancestor_name)
                ancestor = pathseq_df.loc[pathseq_df["name"] == ancestor_name]
            if ancestor.empty:
                ancestor_name = x["taxonomy"].split("|")[-5]
                #print(ancestor_name)
                ancestor = pathseq_df.loc[pathseq_df["name"] == ancestor_name]
        # sometimes multiple levels of taxa have the same name
        if ancestor.shape[0] > 1:
            ancestor = ancestor.loc[ancestor["type"] == parent_taxa[x["type"]]]
        #print(ancestor)
        #print(ancestor.shape)
        if(ancestor.shape[0] == 0):
            return np.nan
        ancestor = ancestor.squeeze()
        #print(x["type"])

        #print(ancestor["type"])
        #assert(ancestor["type"] == parent_taxa[x["type"]])
        return ancestor["tax_id"]


cells = pd.read_csv(snakemake.input[0], sep="\t", dtype={"patient": "str", "sample": "str"})
cells.set_index("cell", inplace=True)
cells = cells.loc[cells.patient == snakemake.wildcards["patient"]]
#tax_level = snakemake.wildcards["tax_level"]
kingdom = snakemake.wildcards["kingdom"]
output = []
for _, cell in cells.iterrows():
    try:
        filename = snakemake.params[0].format(cell.patient, cell["sample"], cell.plate, cell.name)
        pathseq_df = pd.read_csv(filename, sep="\t")
        pathseq_df["sample"] = cell.name
        # for some reason, some species are called "species_group"
        # pathseq_df = pathseq_df.replace({"type": {"species_group": "species"}})
        pathseq_df = pathseq_df.loc[pathseq_df["type"].isin(["superkingdom", "phylum", "class", "order", "family", "genus", "species"])]
        if kingdom == "Viruses":
            pathseq_df = pathseq_df.loc[pathseq_df["taxonomy"].str.startswith("root|Viruses")]
        else:
            pathseq_df = pathseq_df.loc[pathseq_df["kingdom"] == kingdom]
        # pathseq_df = pathseq_df.loc[pathseq_df.unambiguous > 0]
        pathseq_df["parent"] = pathseq_df.apply(get_parent, axis=1)
        edgelist = pathseq_df[["parent", "tax_id"]]
        edgelist = edgelist.rename(columns={"tax_id": "child"})
        output.append(edgelist)

    except OSError as err:
        print("OS error: {0}".format(err))
        print("missing {}".format(cell.name))
        cells.drop(cell.name, inplace=True)
    except:
        print("missing {}".format(cell.name))
        print("Unexpected error:", sys.exc_info()[0])
        raise

df = pd.concat(output)
df = df.drop_duplicates()
df = df.loc[~df.isna().any(axis=1)]
df = df.astype({"parent": int})
df = df.astype({"parent": "str", "child": "str"})
#df["parent"] = df["parent"].apply(lambda x: x.replace("\'", "").replace("#", ""))
#df["child"] = df["child"].apply(lambda x: x.replace("\'", "").replace("#", ""))
df.to_csv(snakemake.output[0], sep="\t", index=False)
