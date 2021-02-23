import pandas as pd
from os.path import join


f = join("data", "{patient}", "{sample}-{plate}-{cell}", "{genome}-paired-count.gff")

col_names = ["seqname", "source", "type", "start", "end", "score", "strand", "frame", "attributes", "read_count"]

cells = pd.read_csv(snakemake.input[0], sep="\t")
genome = snakemake.wildcards["genome"]

rRNA_output = []
protein_coding_output = []

for _, cell in cells.iterrows():
    cell_file = f.format(patient=cell["patient"], sample=cell["sample"], plate=cell["plate"], cell=cell["cell"], genome=genome)
    df = pd.read_csv(cell_file, sep="\t", names=col_names, comment="#")
    df = df.loc[df.type == "gene"]
    df["geneID"] = df["attributes"].apply(lambda x: x.split(";")[0].replace("ID=", ""))
    df = df.set_index(keys="geneID")
    rRNA_df = df.loc[df.attributes.str.contains("gene_biotype=rRNA")]
    rRNA_output.append(rRNA_df["read_count"].rename(cell["cell"]))
    protein_df = df.loc[df.attributes.str.contains("gene_biotype=protein_coding")]
    protein_coding_output.append(protein_df["read_count"].rename(cell["cell"]))

rRNA_df = pd.concat(rRNA_output, axis=1)
rRNA_df.to_csv(snakemake.output[0], sep="\t")
rRNA_df.sum().to_csv(snakemake.output[1], sep="\t", header=None)

protein_coding_df = pd.concat(protein_coding_output, axis=1)
protein_coding_df.to_csv(snakemake.output[2], sep="\t")
protein_coding_df.sum().to_csv(snakemake.output[3], sep="\t", header=None)
