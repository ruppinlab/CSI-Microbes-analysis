
include: "Snakefile"

# Supplementary Tables
TABLE_S1_3 = join("output", "{genome}_gene_reads_supplementary_table1.tsv")
# Files for SRPRISM
PAIRED_SALMONELLA_GFF = join("raw", "SRPRISM", "{patient}", "{sample}-{plate}-{cell}", "{genome}-paired-count.gff")
Salmonella_rRNA_Table = join("output", "Pt0", "S0", "{genome}-cell-by-rRNA-count-table.tsv")
Salmonella_rRNA_read_count = join("output", "Pt0", "S0", "{genome}-rRNA-read-count.tsv")
Salmonella_Gene_Table = join("output", "Pt0", "S0", "{genome}-cell-by-gene-count-table.tsv")
Salmonella_Gene_read_count = join("output", "Pt0", "S0", "{genome}-gene-read-count.tsv")

rule:
    input:
        FIGURE_S1A,
        FIGURE_S1B,
        FIGURE_S1C,
        expand(TABLE_S1_3, genome=["D23580", "LT2"]),

### Figures ###

rule plot_figure_S1:
    input:
        FIGURE_S1A,
        FIGURE_S1B,
        FIGURE_S1C,

### Tables ###

rule generate_Supplementary_Table1_gene_list:
    input:
        "data/units.tsv",
        Salmonella_Gene_Table
    output:
        TABLE_S1_3
    script:
        "src/generate_STable1_gene_read_count.py"

### rules for processing the SRPRISM Salmonella files ###

rule generate_gene_matrices_and_total_counts:
    input:
        "data/units.tsv",
        expand(PAIRED_SALMONELLA_GFF, zip, patient=units["patient"], sample=units["sample"],
               plate=units["plate"], cell=units["cell"], genome=["{genome}"]*units.shape[0])
    output:
        Salmonella_rRNA_Table,
        Salmonella_rRNA_read_count,
        Salmonella_Gene_Table,
        Salmonella_Gene_read_count
    script:
        "src/generate_Salmonella_gene_table.py"


### rules for plotting Supplementary Figure 1 ###

rule plot_figure_S1A:
    input:
        join("data", "units.tsv"),
        Salmonella_Gene_read_count.format(genome="D23580"),
        Salmonella_rRNA_read_count.format(genome="D23580"),
    output:
        FIGURE_S1A
    script:
        "src/plot_figure_S1A.R"

rule plot_figure_S1B:
    input:
        Salmonella_rRNA_read_count.format(genome="D23580"),
        join("data", "units.tsv")
    output:
        FIGURE_S1B
    script:
        "src/plot_figure_S1B.R"

rule plot_figure_S1C:
    input:
        Salmonella_Gene_read_count.format(genome="D23580"),
        join("data", "units.tsv")
    output:
        FIGURE_S1C
    script:
        "src/plot_figure_S1C.R"
