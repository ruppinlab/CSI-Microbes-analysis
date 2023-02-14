import pandas as pd
from numpy import log2


# def get_expected_number_of_double_infected_cells(df):
#     n_infected_by_microbe_df = df.groupby(["patient", "sample"]).apply(lambda x: x.loc[(x[read_df.columns] >= 2).any(axis=1)].shape[0])
#     cells_per_sample = df.groupby(["patient", "sample"]).apply(lambda x: x.shape[0])
#     return (((n_infected_by_microbe_df/cells_per_sample)*(n_infected_by_microbe_df/cells_per_sample))*cells_per_sample).sum()
#
# def get_actual_number_of_double_infected_cells(df):
#     n_infected_by_microbe_df = df.groupby(["patient", "sample"]).apply(lambda x: x.loc[(x[read_df.columns] >= 2).sum(axis=1) >= 2].shape[0])
#     return n_infected_by_microbe_df.sum()
#
# def get_actual_number_of_triple_infected_cells(df):
#     n_infected_by_microbe_df = df.groupby(["patient", "sample"]).apply(lambda x: x.loc[(x[read_df.columns] >= 2).sum(axis=1) >= 3].shape[0])
#     return n_infected_by_microbe_df.sum()

# calculate the
def get_expected_double_infected_cells_per_microbe(df, microbe):
    # how many cells are infected by the microbe of interest?
    n_infected_by_microbe_df = df.groupby(["patient", "sample"]).apply(lambda x: x.loc[(x[microbe] >= 2)].shape[0])
    # how many cells are infected by other microbes?
    other_microbes = list(set(read_df.columns).difference(set([microbe])))
    n_infected_by_other_microbes_df = df.groupby(["patient", "sample"]).apply(lambda x: x.loc[(x[other_microbes] >= 2).any(axis=1)].shape[0])
    cells_per_sample = df.groupby(["patient", "sample"]).apply(lambda x: x.shape[0])
    return (((n_infected_by_microbe_df/cells_per_sample)*(n_infected_by_other_microbes_df/cells_per_sample))*cells_per_sample).sum()

def get_actual_double_infected_cells_per_microbe(df, microbe):
    other_microbes = list(set(read_df.columns).difference(set([microbe])))
    n_infection_df = df.groupby(["patient", "sample"]).apply(lambda x: x.loc[(x[microbe] >= 2) & ((x[other_microbes] >= 2).any(axis=1))].shape[0])
    # print(n_infection_df)
    return n_infection_df.sum()

def get_expected_double_infected_cells_two_microbes(df, microbe1, microbe2):
    # how many cells are infected by the microbe1?
    n_infected_by_microbe1_df = df.groupby(["patient", "sample"]).apply(lambda x: x.loc[(x[microbe1] >= 2)].shape[0])
    # how many cells are infected by microbe2
    n_infected_by_microbe2_df = df.groupby(["patient", "sample"]).apply(lambda x: x.loc[(x[microbe2] >= 2)].shape[0])
    cells_per_sample = df.groupby(["patient", "sample"]).apply(lambda x: x.shape[0])
    return (((n_infected_by_microbe1_df/cells_per_sample)*(n_infected_by_microbe2_df/cells_per_sample))*cells_per_sample).sum()

def get_actual_double_infected_cells_two_microbes(df, microbe1, microbe2):
    n_infection_df = df.groupby(["patient", "sample"]).apply(lambda x: x.loc[(x[microbe1] >= 2) & (x[microbe2] >= 2)].shape[0])
    # print(n_infection_df)
    return n_infection_df.sum()



read_df = pd.read_csv("data/cohort_genera_matrix.tsv", sep="\t", index_col=0)
meta_df = pd.read_csv("data/cohort_metadata.tsv", sep="\t", index_col=0)
df = read_df.merge(meta_df, left_index=True, right_index=True)

output = []

for genera in read_df.columns:
    expected_double_infected_cells_per_microbe = get_expected_double_infected_cells_per_microbe(df, genera)
    double_infected_cells_per_microbe = get_actual_double_infected_cells_per_microbe(df, genera)
    print(genera)
    print(expected_double_infected_cells_per_microbe)
    print(double_infected_cells_per_microbe)
    print(log2(double_infected_cells_per_microbe/expected_double_infected_cells_per_microbe))
    output.append(pd.Series(data={"microbe": genera,
                                  "expected_n_double_infected": expected_double_infected_cells_per_microbe,
                                  "actual_n_double_infected": double_infected_cells_per_microbe,
                                  "log2FC": log2(double_infected_cells_per_microbe/expected_double_infected_cells_per_microbe)
                                  }))

output_df = pd.concat(output, axis=1).T
output_df.sort_values(by="log2FC")

output = []

for genera1 in read_df.columns:
    for genera2 in read_df.columns:
        if genera1 == genera2:
            continue
        expected_double_infected_cells_per_microbe = get_expected_double_infected_cells_two_microbes(df, genera1, genera2)
        double_infected_cells_per_microbe = get_actual_double_infected_cells_two_microbes(df, genera1, genera2)
        output.append(pd.Series(data={"microbe1": genera1,
                                      "microbe2": genera2,
                                      "expected_n_double_infected": expected_double_infected_cells_per_microbe,
                                      "actual_n_double_infected": double_infected_cells_per_microbe,
                                      "log2FC": log2(double_infected_cells_per_microbe/expected_double_infected_cells_per_microbe)
                                      }))


output_df = pd.concat(output, axis=1).T
output_df.sort_values(by="log2FC")
output_df.loc[output_df.notna().all(axis=1)].sort_values(by="actual_n_double_infected")

myeloid_df = df.loc[df.celltype1 == "Myeloid"]
output = []

for genera in read_df.columns:
    expected_double_infected_cells_per_microbe = get_expected_double_infected_cells_per_microbe(myeloid_df, genera)
    double_infected_cells_per_microbe = get_actual_double_infected_cells_per_microbe(myeloid_df, genera)
    print(log2(double_infected_cells_per_microbe/expected_double_infected_cells_per_microbe))
    output.append(pd.Series(data={"microbe": genera,
                                  "expected_n_double_infected": expected_double_infected_cells_per_microbe,
                                  "actual_n_double_infected": double_infected_cells_per_microbe,
                                  "log2FC": log2(double_infected_cells_per_microbe/expected_double_infected_cells_per_microbe)
                                  }))

output_df = pd.concat(output, axis=1).T
output_df.sort_values(by="log2FC")
# non_myeloid_df = df.loc[df.celltype1 != "Myeloid"]
#
# for genera in read_df.columns:
#     expected_double_infected_cells_per_celltype = get_expected_infected_cells(non_myeloid_df, genera)
#     infected_cells_per_celltype = get_n_infected_cells(non_myeloid_df, genera)
#     print(genera)
#     print(expected_infected_cells_per_celltype)
#     print(log2(infected_cells_per_celltype/expected_infected_cells_per_celltype))
