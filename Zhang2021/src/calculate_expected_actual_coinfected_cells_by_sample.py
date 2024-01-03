import pandas as pd
import numpy as np
import math
from scipy.stats import hypergeom, combine_pvalues
from statsmodels.stats.multitest import fdrcorrection


def get_expected_coinfected_cells(df, genera1, other_genera, cutoff):
    # calculate the expected number of infected cells per cell-type
    n_expected_df = df.groupby(["patient", "sample"]).apply(lambda x: ((x.loc[(x[genera1] >= cutoff)].shape[0]/x.shape[0])*(x.loc[(x[other_genera] >= cutoff).any(axis=1)].shape[0]/x.shape[0]))*x.shape[0])
    n_expected_df.name = "n_expected_infected"
    n_expected_df = n_expected_df.reset_index()
    return n_expected_df


# def get_n_coinfected_cells(df, genera1, genera2, cutoff):
#     # calculate the actual number of infected cells per cell-type
#     n_infection_df = df.groupby(["patient", "sample"]).apply(lambda x: x.loc[(x[genera1] >= cutoff) & (x[genera2] >= cutoff)].shape[0])
#     return n_infection_df

# first step is to calculate a p-value for each sample for each cell-type
# we were previously calculating this for each genera in each cell-type
# let's calculate this using the hypergeometric enrichment test
def single_sample_g1_g2_hypergeom(df, sample, g1, other_genera, cutoff):
    sample_df = df.loc[df["sample"] == sample]
    #celltype_of_interest_df = sample_df.loc[sample_df[celltype_col] == celltype_of_interest]
    #non_celltype_of_interest_df = sample_df.loc[sample_df[celltype_col] != celltype_of_interest]
    # get the number of infected cells from celltype of interest
    n = sample_df.loc[sample_df[g1] >= cutoff].shape[0]
    N = sample_df.loc[(sample_df[other_genera] >= cutoff).any(axis=1)].shape[0]
    k = sample_df.loc[(sample_df[g1] >= cutoff) & (sample_df[other_genera] >= cutoff).any(axis=1)].shape[0]
    M = sample_df.shape[0]
    hpd = hypergeom(M, n, N)
    p = hpd.pmf(k)
    return p


def calculate_p_value(df, g1, other_genera, cutoff):
    samples = df["sample"].unique()
    output = []
    for sample in samples:
        pval = single_sample_g1_g2_hypergeom(df, sample, g1, other_genera, cutoff)
        output.append(pd.Series({"sample": sample, "g1": g1, "pvalue": pval}))
    return pd.concat(output, axis=1).T


def calculate_coocurrence_pvalues(df, g1, other_genera, cutoff):
    #expected_coinfected_cells = get_expected_coinfected_cells(df, g1, other_genera, cutoff)
    coinfection_pvalues = calculate_p_value(df, g1, other_genera, cutoff)
    #coinfection_pvalue_df = expected_coinfected_cells[["sample", "n_expected_infected"]].merge(coinfection_pvalues[["sample", "pvalue"]], on="sample")
    return coinfection_pvalues
    # Stouffer's returns -inf if any p-value is equal to 1 so we need to set .9999 as the maximum value
    # coinfection_pvalue_df["pvalue"] = coinfection_pvalue_df["pvalue"].apply(lambda x: min(x, .9999))
    # return combine_pvalues(coinfection_pvalue_df["pvalue"].astype("float").values, method="stouffer", weights=coinfection_pvalue_df["n_expected_infected"].astype("float").values)[1]



# units = pd.read_csv("data/units.tsv", sep="\t")
# sample_df = units[["patient", "sample"]].drop_duplicates()
# sample_df = sample_df.loc[sample_df["sample"].str.contains("T")]
# sample_df = sample_df.loc[~sample_df.patient.isin(["P38", "P39", "P40"])]
# from os.path import join
# SAMPLE_MICROBE_READ_TABLE = join("output", "{patient}", "{sample}", "{tax_level}_{method}_{kingdom}_reads.tsv")
# files = [SAMPLE_MICROBE_READ_TABLE.format(patient=row.patient, sample=row["sample"], tax_level="genus", method="PathSeq", kingdom="All") for _, row in sample_df.iterrows()]
files = snakemake.input["microbe_reads"]
output = []
for f in files:
    df = pd.read_csv(f, sep="\t", index_col=0)
    output.append(df)

# concat merges on the index by default
read_df = pd.concat(output, join="outer", axis=1).fillna(0)

# PATIENT_PATHSEQ_TAXID_MAP = join("output", "{patient}", "tax_id_map_{kingdom}_{method}.tsv")
# patient_df = units[["patient"]].drop_duplicates()
# patient_df = patient_df.loc[~patient_df.patient.isin(["P38", "P39", "P40"])]
# tax_files = [PATIENT_PATHSEQ_TAXID_MAP.format(patient=row.patient, tax_level="genus", method="PathSeq", kingdom="All") for _, row in sample_df.iterrows()]
# add the tax id information
tax_files = snakemake.input["tax_ids"]
output = []
for f in tax_files:
    output.append(pd.read_csv(f, sep="\t"))

tax_map = pd.concat(output).drop_duplicates()

read_df = read_df.merge(tax_map, left_index=True, right_on="tax_id")
read_df = read_df.set_index("name")
read_df = read_df.drop(columns=["tax_id", "taxa_level"])
read_df = read_df.T
genera = read_df.columns

meta_df = pd.read_csv("data/units.tsv", sep="\t")
meta_df = meta_df.loc[meta_df["sample"].str.contains("T")]
meta_df.index = meta_df.apply(lambda x: "{}-{}".format(x["sample"], x["barcode"]), axis=1)
df = meta_df[["patient", "sample", "barcode", "celltype1", "celltype2"]].merge(read_df, left_index=True, right_index=True)

cutoff = int(snakemake.wildcards["cutoff"])

# just focus on genera with >= 2 UMIs in any cell
infecting_genera = genera[(df[genera] >= cutoff).sum() > 10]


# now let's calculate a co-occurrence row for all bacteria
# key thing is that we need to correctly calculate the expected co-infected
# co-infected cells is jut the number of cells with >= 2 genera above the treshold
output = []
comparison_genera = infecting_genera.tolist()
for g1 in infecting_genera:
    comparison_genera.remove(g1)
    expected_coinfected_cells = get_expected_coinfected_cells(df, g1, comparison_genera, cutoff)
    # print(expected_coinfected_cells)
    if expected_coinfected_cells["n_expected_infected"].sum() == 0:
        continue
    else:
        pval = calculate_coocurrence_pvalues(df, g1, comparison_genera, cutoff)
        expected_coinfected_cells = expected_coinfected_cells.merge(pval, on=["sample"])
    output.append(expected_coinfected_cells)

co_df = pd.concat(output)

co_df["pvalue"] = co_df["pvalue"].apply(lambda x: min(x, .9999))
pvalue_by_sample = co_df.groupby("sample").apply(lambda x: combine_pvalues(x["pvalue"].astype("float"), method="stouffer", weights=x["n_expected_infected"].astype("float").values)[1])
expected_coinfection_by_sample = co_df.groupby("sample")["n_expected_infected"].sum()
pvalue_by_sample
pvalue = combine_pvalues(pvalue_by_sample.fillna(.9999), method="stouffer", weights=expected_coinfection_by_sample)[1]
# output.append(pd.Series({"g1": g1, "g2": g2, "expected_co_infected": expected_coinfected_cells["n_expected_infected"].sum(), "co_infected": coinfected_cells.sum(), "pvalue": pval}))


sample_expected_coinfected = expected_coinfection_by_sample.sum()
sample_coinfected = df.groupby(["patient", "sample"]).apply(lambda x: ((x[infecting_genera] >= cutoff).sum(axis=1) >= 2).sum())



pd.DataFrame({"coinfected": [sample_coinfected.sum()], "expected_coinfected": [sample_expected_coinfected], "log2_FC": [math.log2(sample_coinfected.sum()/sample_expected_coinfected)], "pvalue": [pvalue]}).to_csv(snakemake.output[0], index=False, sep="\t")
