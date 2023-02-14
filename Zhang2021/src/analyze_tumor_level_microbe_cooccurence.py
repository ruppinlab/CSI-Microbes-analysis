import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import fdrcorrection


df = pd.read_csv("data/units.tsv", sep="\t")
# focus only on tumor samples
df = df.loc[df["sample"].str.contains("T")]

tax_level = "genus"
read_file = "output/{}/{}_PathSeq_All_reads.tsv"
output = []
for p in df.patient.unique():
    output.append(pd.read_csv(read_file.format(p, tax_level), sep="\t", index_col=0))

read_df = pd.concat(output, join="outer", axis=1).fillna(0)
read_df = read_df.T

meta_df = pd.read_csv("data/units.tsv", sep="\t")
meta_df.index = meta_df.apply(lambda x: "{}-{}".format(x["sample"], x["barcode"]), axis=1)

tax_map_file = "output/{}/tax_id_map_All_PathSeq.tsv"
output = []
for p in df.patient.unique():
    df = pd.read_csv(tax_map_file.format(p), sep="\t", index_col=0)
    output.append(df)

tax_df = pd.concat(output).drop_duplicates()
tax_df = tax_df.loc[tax_df["taxa_level"] == tax_level]
d = dict(zip(tax_df["tax_id"], tax_df.index))
read_df = read_df.rename(columns=d)
df = meta_df.merge(read_df, left_index=True, right_index=True)

df = df.loc[df["sample"].str.contains("T")]

min_umis = 2
genera = read_df.columns[(df[read_df.columns] >= min_umis).sum() >= 10]
patient_genera_matrix = df.groupby(["patient"]).apply(lambda x: (x[genera] >= min_umis).sum())
n_patients = patient_genera_matrix.shape[0]
min_cells = 1
output = []
for g1 in genera:
    for g2 in genera:
        # number of patients where g1 is found
        g1_n_patients = patient_genera_matrix.loc[patient_genera_matrix[g1] >= min_cells].shape[0]
        # number of patients where g2 is found
        g2_n_patients = patient_genera_matrix.loc[patient_genera_matrix[g2] >= min_cells].shape[0]
        overlap_n_patients = patient_genera_matrix.loc[(patient_genera_matrix[g1] >= min_cells) & (patient_genera_matrix[g2] >= min_cells)].shape[0]
        #print("co-occurrence p-value for {} and {}".format(b1, b2))
        pval = hypergeom.sf(overlap_n_patients-1, n_patients, g1_n_patients, g2_n_patients)
        output.append(pd.Series(data={"g1": g1, "g2": g2, "pval": pval, "g1_n": g1_n_patients, "g2_n": g2_n_patients, "overlap_n_patients": overlap_n_patients}))


overlap_df = pd.concat(output, axis=1).T
overlap_df = overlap_df.loc[overlap_df["g1"] != overlap_df["g2"]]
overlap_df = overlap_df.loc[overlap_df["g2"] > overlap_df["g1"]]
overlap_df["fdr"] = fdrcorrection(overlap_df["pval"], method="i")[1]
overlap_df = overlap_df.sort_values(by="pval")
