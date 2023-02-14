
include: "Snakefile"

# Supplementary Tables
TABLE_S1_3 = join("output", "{genome}_gene_reads_supplementary_table1.tsv")
# Files for SRPRISM
PAIRED_GFF = join("raw", "SRPRISM", "{patient}", "{sample}-{plate}-{cell}", "{genome}-paired-count.gff")
read_type_table = join("output", "Pt0", "S0", "{genome}-read-type-table.tsv")
# Salmonella_rRNA_read_count = join("output", "Pt0", "S0", "{genome}-rRNA-read-count.tsv")
# Salmonella_Gene_Table = join("output", "Pt0", "S0", "{genome}-cell-by-gene-count-table.tsv")
# Salmonella_Gene_read_count = join("output", "Pt0", "S0", "{genome}-gene-read-count.tsv")

rule:
    input:
        expand("output/plots/piechart_RNA_type_{genome}.svg", genome="Fn")

### Figures ###

# rule plot_figure_S1:
#     input:
#         FIGURE_S1A,
#         FIGURE_S1B,
#         FIGURE_S1C,

### Tables ###

rule plot_piechart_RNA_class:
    input:
        read_type_table
    output:
        "output/plots/piechart_RNA_type_{genome}.svg"
    script:
        "src/plot_piechart_RNA_class.py"

### rules for processing the SRPRISM Salmonella files ###
rule generate_gene_matrices_and_total_counts:
    input:
        "data/units.tsv",
        expand(PAIRED_GFF, zip, patient=units["patient"], sample=units["sample"],
               plate=units["plate"], cell=units["cell"], genome=["{genome}"]*units.shape[0])
    output:
        read_type_table
    script:
        "src/generate_Fusobacterium_gene_tables_from_SRPRISM_SS2.py"


# rsync -avc --include='*-paired-count.gff' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Robinson2022-SS2/output/ raw/
# rsync -avc --include='Fn_read_counts.tsv' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Robinson2022-SS2/output/ raw/
