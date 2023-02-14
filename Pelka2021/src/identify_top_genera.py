import pandas as pd
from scipy.stats import hypergeom


# loop through the superkingdom reads
# 2 = "Bacteria"
# 2759 = "Eukaryota"
# 10239 = "Viruses"
files = snakemake.input["microbe_reads"]
output = []
for f in files:
    df = pd.read_csv(f, sep="\t", index_col=0)
    output.append(df)
# concat merges on the index by default
read_df = pd.concat(output, join="outer", axis=1).fillna(0)
read_df = read_df.T

# get the tax_dfs
files = snakemake.input["taxid_maps"]
output = []
for f in files:
    df = pd.read_csv(f, sep="\t", index_col=0)
    output.append(df)

tax_df = pd.concat(output).drop_duplicates()
tax_df = tax_df.loc[tax_df["taxa_level"] == snakemake.wildcards["tax_level"]]
d = dict(zip(tax_df["tax_id"], tax_df.index))
read_df = read_df.rename(columns=d)

meta_df = pd.read_csv(snakemake.input["meta_data"], sep="\t")
meta_df.index = meta_df.apply(lambda x: "{}-{}".format(x["sample"], x["barcode"]), axis=1)
df = read_df.merge(meta_df, left_index=True, right_index=True)

# only focus on cells from the tumor
df = df.loc[df["Is_Tumor"] == "Yes"]
genera_with_min_umis = read_df.columns[(df[read_df.columns] >= 2).any()]

columns = list(meta_df.columns) + list(genera_with_min_umis)
df = df[columns]

output = []
for taxa in list(genera_with_min_umis):
    output.append(pd.DataFrame(data={"taxa": [taxa],
                       "n_umis": [df[taxa].sum()],
                       "n_cells_min_1_umi": [df.loc[df[taxa] >= 1].shape[0]],
                       "n_cells_min_2_umi": [df.loc[df[taxa] >= 2].shape[0]],
                       "n_patients_min_1_umi": [df.loc[df[taxa] >= 1].patient.nunique()],
                       "n_patients_min_2_umi": [df.loc[df[taxa] >= 2].patient.nunique()]
                       }))

genera_df = pd.concat(output)
# filter out mycoplasma
genera_df = genera_df.loc[genera_df["taxa"] != "Mycoplasma"]
# focus on genera with at least 10 "infected cells"
# genera_df = genera_df.loc[genera_df["n_cells_min_2_umi"] > 10]
#
genera_df.sort_values(by="n_umis", ascending=False).to_csv(snakemake.output[0], sep="\t", index=False)


# min_umis = 2
# top_bacteria = genera_df["taxa"]
# output = []
# for b in top_bacteria:
#     for p in df.patient.unique():
#         # patient_df = df.loc[(df.patient == p) & (df["celltype1"] != "Epi")]
#         patient_df = df.loc[df.patient == p]
#         if (patient_df[b] >= min_umis).any():
#             myeloid_cells = patient_df.loc[patient_df["celltype1"].isin(["Myeloid"])]
#             myeloid_pos = myeloid_cells.loc[myeloid_cells[b] >= min_umis]
#             pos_cells = patient_df.loc[patient_df[b] >= min_umis]
#             pval = hypergeom.sf(myeloid_pos.shape[0]-1, patient_df.shape[0], myeloid_cells.shape[0], pos_cells.shape[0])
#             output.append(pd.Series(data={"bacteria": b,
#                                           "patient": p,
#                                           "pval": pval,
#                                           "pos_cells": pos_cells.shape[0],
#                                           "myeloid_pos_cells": myeloid_pos.shape[0],
#                                           "myeloid_cells": myeloid_cells.shape[0],
#                                           "total_cells": patient_df.shape[0],
#                                           "min_umis": min_umis,
#                                           }))
#
# mye_df = pd.concat(output, axis=1).T.sort_values(by="pval")
#
# min_umis = 1
# top_bacteria = genera_df["taxa"]
# output = []
# for b in top_bacteria:
#     for p in df.patient.unique():
#         # patient_df = df.loc[(df.patient == p) & (df["celltype1"] != "Epi")]
#         patient_df = df.loc[df.patient == p]
#         if (patient_df[b] >= min_umis).any():
#             myeloid_cells = patient_df.loc[patient_df["celltype1"].isin(["Myeloid", "Epi"])]
#             myeloid_pos = myeloid_cells.loc[myeloid_cells[b] >= min_umis]
#             pos_cells = patient_df.loc[patient_df[b] >= min_umis]
#             pval = hypergeom.sf(myeloid_pos.shape[0]-1, patient_df.shape[0], myeloid_cells.shape[0], pos_cells.shape[0])
#             output.append(pd.Series(data={"bacteria": b,
#                                           "patient": p,
#                                           "pval": pval,
#                                           "pos_cells": pos_cells.shape[0],
#                                           "myeloid_pos_cells": myeloid_pos.shape[0],
#                                           "myeloid_cells": myeloid_cells.shape[0],
#                                           "total_cells": patient_df.shape[0],
#                                           "min_umis": min_umis,
#                                           }))
#
# mye_epi_df = pd.concat(output, axis=1).T.sort_values(by="pval")
#
#
# min_umis = 1
# top_bacteria = genera_df["taxa"]
# output = []
# for b in top_bacteria:
#     for p in df.patient.unique():
#         patient_df = df.loc[(df.patient == p) & (df["celltype1"] != "Myeloid")]
#         # patient_df = df.loc[(df.patient == p)]
#         if (patient_df[b] >= min_umis).any():
#             myeloid_cells = patient_df.loc[patient_df["celltype1"] == "Epi"]
#             myeloid_pos = myeloid_cells.loc[myeloid_cells[b] >= min_umis]
#             pos_cells = patient_df.loc[patient_df[b] >= min_umis]
#             pval = hypergeom.sf(myeloid_pos.shape[0]-1, patient_df.shape[0], myeloid_cells.shape[0], pos_cells.shape[0])
#             output.append(pd.Series(data={"bacteria": b,
#                                           "patient": p,
#                                           "pval": pval,
#                                           "pos_cells": pos_cells.shape[0],
#                                           "myeloid_pos_cells": myeloid_pos.shape[0],
#                                           "myeloid_cells": myeloid_cells.shape[0],
#                                           "total_cells": patient_df.shape[0],
#                                           "min_umis": min_umis,
#                                           }))
#
# pd.concat(output, axis=1).T.sort_values(by="pval")
#
#
# min_umis = 1
# top_bacteria = genera_df["taxa"]
# output = []
# for b in top_bacteria:
#     patients = df.loc[df[b] >= min_umis]["patient"].unique()
#     patient_df = df.loc[df.patient.isin(patients)] # == p) & (df["celltype1"] != "Myeloid")]
#     myeloid_cells = patient_df.loc[patient_df["celltype1"] == "Myeloid"]
#     myeloid_pos = myeloid_cells.loc[myeloid_cells[b] >= min_umis]
#     pos_cells = patient_df.loc[patient_df[b] >= min_umis]
#     pval = hypergeom.sf(myeloid_pos.shape[0]-1, patient_df.shape[0], myeloid_cells.shape[0], pos_cells.shape[0])
#     output.append(pd.Series(data={"bacteria": b,
#                                   "patient": p,
#                                   "pval": pval,
#                                   "pos_cells": pos_cells.shape[0],
#                                   "myeloid_pos_cells": myeloid_pos.shape[0],
#                                   "myeloid_cells": myeloid_cells.shape[0],
#                                   "total_cells": patient_df.shape[0],
#                                   "min_umis": min_umis,
#                                   }))
#
# pd.concat(output, axis=1).T.sort_values(by="pval")
#
# min_umis = 1
# top_bacteria = genera_df["taxa"]
# output = []
# for b in top_bacteria:
#     patients = df.loc[df[b] >= min_umis]["patient"].unique()
#     patient_df = df.loc[df.patient.isin(patients)] # == p) & (df["celltype1"] != "Myeloid")]
#     myeloid_cells = patient_df.loc[patient_df["celltype1"].isin(["B", "Plasma"])]
#     myeloid_pos = myeloid_cells.loc[myeloid_cells[b] >= min_umis]
#     pos_cells = patient_df.loc[patient_df[b] >= min_umis]
#     pval = hypergeom.sf(myeloid_pos.shape[0]-1, patient_df.shape[0], myeloid_cells.shape[0], pos_cells.shape[0])
#     output.append(pd.Series(data={"bacteria": b,
#                                   "patient": p,
#                                   "pval": pval,
#                                   "pos_cells": pos_cells.shape[0],
#                                   "myeloid_pos_cells": myeloid_pos.shape[0],
#                                   "myeloid_cells": myeloid_cells.shape[0],
#                                   "total_cells": patient_df.shape[0],
#                                   "min_umis": min_umis,
#                                   }))
#
# pd.concat(output, axis=1).T.sort_values(by="pval")
#
# min_umis = 2
# top_bacteria = genera_df["taxa"]
# top_bacteria = ["Lymphocryptovirus"]
# output = []
# for b in top_bacteria:
#     for p in df.patient.unique():
#         print(p)
#         patient_df = df.loc[df.patient == p]
#         print((patient_df[b] >= min_umis).any())
#         if (patient_df[b] >= min_umis).any():
#             myeloid_cells = patient_df.loc[patient_df["celltype1"].isin(["B", "Plasma"])]
#             myeloid_pos = myeloid_cells.loc[myeloid_cells[b] >= min_umis]
#             pos_cells = patient_df.loc[patient_df[b] >= min_umis]
#             pval = hypergeom.sf(myeloid_pos.shape[0]-1, patient_df.shape[0], myeloid_cells.shape[0], pos_cells.shape[0])
#             output.append(pd.Series(data={"bacteria": b,
#                                           "patient": p,
#                                           "pval": pval,
#                                           "pos_cells": pos_cells.shape[0],
#                                           "myeloid_pos_cells": myeloid_pos.shape[0],
#                                           "myeloid_cells": myeloid_cells.shape[0],
#                                           "total_cells": patient_df.shape[0],
#                                           "min_umis": min_umis,
#                                           }))
#
# pd.concat(output, axis=1).T.sort_values(by="pval")
#
#
#         b1_n = df.loc[df[b1] >= min_umis].shape[0]
#         b2_n = df.loc[df[b2] >= min_umis].shape[0]
#         overlap_n = df.loc[(df[b1] >= min_umis) & (df[b2] >= min_umis)].shape[0]
#         b1_n_patients = df.loc[df[b1] >= min_umis]["patient"].nunique()
#         b2_n_patients = df.loc[df[b2] >= min_umis]["patient"].nunique()
#         n_patients_overlap = df.loc[(df[b1] >= min_umis) & (df[b2] >= min_umis)]["patient"].nunique()
#         pval = hypergeom.sf(overlap_n-1, df.shape[0], b1_n, b2_n)
#         output.append(pd.Series(data={"b1": b1,
#                                       "b2": b2,
#                                       "pval": pval,
#                                       "b1_n": b1_n,
#                                       "b2_n": b2_n,
#                                       "overlap_n": overlap_n,
#                                       "b1_n_patients": b1_n_patients,
#                                       "b2_n_patients": b2_n_patients,
#                                       "n_patients_overlap": n_patients_overlap
#                                       }))
#
#
# min_umis = 2
# top_bacteria = genera_df["taxa"]
# output = []
# for b1 in top_bacteria:
#     for b2 in top_bacteria:
#         b1_n = df.loc[df[b1] >= min_umis].shape[0]
#         b2_n = df.loc[df[b2] >= min_umis].shape[0]
#         overlap_n = df.loc[(df[b1] >= min_umis) & (df[b2] >= min_umis)].shape[0]
#         b1_n_patients = df.loc[df[b1] >= min_umis]["patient"].nunique()
#         b2_n_patients = df.loc[df[b2] >= min_umis]["patient"].nunique()
#         n_patients_overlap = df.loc[(df[b1] >= min_umis) & (df[b2] >= min_umis)]["patient"].nunique()
#         pval = hypergeom.sf(overlap_n-1, df.shape[0], b1_n, b2_n)
#         output.append(pd.Series(data={"b1": b1,
#                                       "b2": b2,
#                                       "pval": pval,
#                                       "b1_n": b1_n,
#                                       "b2_n": b2_n,
#                                       "overlap_n": overlap_n,
#                                       "b1_n_patients": b1_n_patients,
#                                       "b2_n_patients": b2_n_patients,
#                                       "n_patients_overlap": n_patients_overlap
#                                       }))
#
# no_163_169_df = df.loc[~df.patient.isin(["C163", "C169"])]
#
# output = []
# min_umis = 1
# for b1 in top_bacteria:
#     for b2 in top_bacteria:
#         b1_n = no_163_169_df.loc[no_163_169_df[b1] >= min_umis].shape[0]
#         b2_n = no_163_169_df.loc[no_163_169_df[b2] >= min_umis].shape[0]
#         overlap_n = no_163_169_df.loc[(no_163_169_df[b1] >= min_umis) & (no_163_169_df[b2] >= min_umis)].shape[0]
#         b1_n_patients = no_163_169_df.loc[no_163_169_df[b1] >= min_umis]["patient"].nunique()
#         b2_n_patients = no_163_169_df.loc[no_163_169_df[b2] >= min_umis]["patient"].nunique()
#         n_patients_overlap = no_163_169_df.loc[(no_163_169_df[b1] >= min_umis) & (no_163_169_df[b2] >= min_umis)]["patient"].nunique()
#         pval = hypergeom.sf(overlap_n-1, df.shape[0], b1_n, b2_n)
#         output.append(pd.Series(data={"b1": b1,
#                                       "b2": b2,
#                                       "pval": pval,
#                                       "b1_n": b1_n,
#                                       "b2_n": b2_n,
#                                       "overlap_n": overlap_n,
#                                       "b1_n_patients": b1_n_patients,
#                                       "b2_n_patients": b2_n_patients,
#                                       "n_patients_overlap": n_patients_overlap
#                                       }))
#
#
# overlap_df = pd.concat(output, axis=1).T
# overlap_df = overlap_df.loc[overlap_df["b1"] != overlap_df["b2"]]
# overlap_df.loc[overlap_df["b2"] > overlap_df["b1"]]
#
# # overlap_df = overlap_df.drop_duplicates(subset=["b1", "b2"])
# overlap_df["fdr"] = fdrcorrection(overlap_df["pval"], method="i")[1]
# overlap_df = overlap_df.sort_values(by="pval")
# print(overlap_df)
