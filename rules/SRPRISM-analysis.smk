
PAIRED_GFF = join("raw", "SRPRISM", "{patient}", "{sample}-{plate}-{cell}", "{genome}-paired-count.gff")
rRNA_table = join("output", "{patient}", "{genome}-cell-by-rRNA-count-table.tsv")
rRNA_read_count = join("output", "{patient}", "{genome}-rRNA-read-count.tsv")
Gene_table = join("output", "{patient}", "{genome}-cell-by-gene-count-table.tsv")
Gene_read_count = join("output", "{patient}", "{genome}-gene-read-count.tsv")

gene_level_table = join("output", "{patient}", "{genome}-gene-read-table.tsv")
type_level_table = join("output", "{patient}", "{genome}-type-read-table.tsv")

# print(Gene_read_count)

def get_input_gff_files(wildcards):
    my_df = units.loc[units.patient == wildcards.patient]
    return expand(PAIRED_GFF, zip, patient=my_df["patient"], sample=my_df["sample"],
           plate=my_df["plate"], cell=my_df["cell"], genome=[wildcards.genome]*my_df.shape[0])


rule generate_Salmonella_gene_matrices_and_total_counts:
    wildcard_constraints:
        genome="LT2|D23580"
    input:
        "data/units.tsv",
        get_input_gff_files
    output:
        rRNA_table,
        rRNA_read_count,
        Gene_table,
        Gene_read_count
    script:
        "../src/generate_Salmonella_gene_tables_from_SRPRISM.py"

rule generate_FN_gene_matrices_and_total_counts:
    wildcard_constraints:
        genome="FN"
    input:
        "data/units.tsv",
        get_input_gff_files
    output:
        gene_level_table,
        type_level_table,
    script:
        "../src/generate_Fusobacterium_gene_tables_from_SRPRISM.py"
