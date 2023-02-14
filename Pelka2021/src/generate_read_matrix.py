import pandas as pd

df = pd.read_csv("data/units.tsv", sep="\t")
df = df.loc[df["Is_Tumor"] == "Yes"]
df = df.loc[df["10x_chemistry"] == "v3"]

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

genera_of_focus = ['Alloprevotella', 'Bacteroides', 'Campylobacter', 'Capnocytophaga',
                   'Corynebacterium', 'Fusobacterium', 'Leptotrichia', 'Mycoplasma',
                   'Parabacteroides', 'Phocaeicola', 'Prevotella', 'Pseudomonas',
                   'Serratia', 'Sphingomonas', 'Streptococcus', 'Treponema']
genera_of_focus = list(set(genera_of_focus).intersection(read_df.columns))

read_df = read_df[genera_of_focus]

df = meta_df.merge(read_df, left_index=True, right_index=True)
df = df.loc[df["Is_Tumor"] == "Yes"]
df = df.loc[df["10x_chemistry"] == "v3"]

df[read_df.columns].to_csv("output/cohort_genera_matrix.tsv", sep="\t")
df[meta_df.columns].to_csv("output/cohort_metadata.tsv", sep="\t")



# tumor_df = df.loc[df.Tumor == "Tumor"]
# tumor_df = tumor_df.loc[tumor_df.celltype2.str.startswith("Epi")]
# tumor_df = tumor_df.loc[tumor_df.celltype2 != "Epithelial"]
#
# top_genera = read_df.columns[tumor_df[read_df.columns].sum() > 50]
# min_umis = 1
# microbe = "Pseudomonas"
# output = []
# for celltype in tumor_df.celltype2.unique():
#     n_microbe_positive = tumor_df.loc[tumor_df[microbe] >= min_umis].shape[0]
#     n_celltype = tumor_df.loc[tumor_df["celltype2"] == celltype].shape[0]
#     n_celltype_microbe_positive = tumor_df.loc[(tumor_df["celltype2"] == celltype) & (tumor_df[microbe] >= min_umis)].shape[0]
#     pval = hypergeom.sf(n_celltype_microbe_positive-1, tumor_df.shape[0], n_microbe_positive, n_celltype)
#     output.append(pd.Series(data={"microbe": microbe, "celltype": celltype,
#                                   "pval": pval,
#                                   "n_microbe_positive": n_microbe_positive,
#                                   "n_celltype": n_celltype, "n_celltype_microbe_positive": n_celltype_microbe_positive}))
#
# min_umis = 1
# top_genera = read_df.columns[tumor_df[read_df.columns].sum() > 50]
# output = []
# for celltype in tumor_df.celltype2.unique():
#     for patient in tumor_df.patient.unique():
#         for microbe in top_genera:
#             patient_df = tumor_df.loc[tumor_df.patient == patient]
#             n_microbe_positive = patient_df.loc[patient_df[microbe] >= min_umis].shape[0]
#             n_celltype = patient_df.loc[patient_df["celltype2"] == celltype].shape[0]
#             n_celltype_microbe_positive = patient_df.loc[(patient_df["celltype2"] == celltype) & (patient_df[microbe] >= min_umis)].shape[0]
#             pval = hypergeom.sf(n_celltype_microbe_positive-1, patient_df.shape[0], n_microbe_positive, n_celltype)
#             output.append(pd.Series(data={"microbe": microbe, "patient": patient,
#                                           "celltype": celltype, "pval": pval,
#                                           "n_microbe_positive": n_microbe_positive,
#                                           "n_celltype": n_celltype,
#                                           "n_celltype_microbe_positive": n_celltype_microbe_positive}))
#
# top_genera = read_df.columns[tumor_df[read_df.columns].sum() > 50]
# min_umis = 1
# output = []
# for celltype in tumor_df.celltype2.unique():
#     for microbe in top_genera:
#         n_microbe_positive = tumor_df.loc[tumor_df[microbe] >= min_umis].shape[0]
#         n_celltype = tumor_df.loc[tumor_df["celltype2"] == celltype].shape[0]
#         n_celltype_microbe_positive = tumor_df.loc[(tumor_df["celltype2"] == celltype) & (tumor_df[microbe] >= min_umis)].shape[0]
#         pval = hypergeom.sf(n_celltype_microbe_positive-1, tumor_df.shape[0], n_microbe_positive, n_celltype)
#         output.append(pd.Series(data={"microbe": microbe, "celltype": celltype,
#                                       "pval": pval,
#                                       "n_microbe_positive": n_microbe_positive,
#                                       "n_celltype": n_celltype, "n_celltype_microbe_positive": n_celltype_microbe_positive}))
#
#
# min_umis = 2
# output = []
# for celltype in tumor_df.celltype2.unique():
#     n_microbe_positive = tumor_df.loc[(tumor_df[read_df.columns] >= min_umis).any(axis=1)].shape[0]
#     n_celltype = tumor_df.loc[tumor_df["celltype2"] == celltype].shape[0]
#     n_celltype_microbe_positive = tumor_df.loc[(tumor_df["celltype2"] == celltype) & (tumor_df[read_df.columns] >= min_umis).any(axis=1)].shape[0]
#     pval = hypergeom.sf(n_celltype_microbe_positive-1, tumor_df.shape[0], n_microbe_positive, n_celltype)
#     output.append(pd.Series(data={"celltype": celltype,
#                                   "pval": pval,
#                                   "n_microbe_positive": n_microbe_positive,
#                                   "n_celltype": n_celltype, "n_celltype_microbe_positive": n_celltype_microbe_positive}))
#
#
#
# for genera in genera_of_interest:
#     print(df.groupby("Smoking status")[genera].mean())
#     print(df.groupby("Smoking status")[genera].median())
#     print(ranksums(df.loc[df["Smoking status"] == "Smoker"][genera], df.loc[df["Smoking status"] == "Non-smoker"][genera]))
