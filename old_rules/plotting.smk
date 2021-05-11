# output plots
# total read box plots
LOG2_TOTAL_READS_BOXPLOT = join("output", "{tax_level}-{method}-{celltype}-log2-total-microbial-reads-boxplot.png")
TOTAL_READS_BOXPLOT = join("output", "{tax_level}-{method}-{celltype}-total-microbial-reads-boxplot.png")
# microbial specific box plots
LOG2_MICROBE_READ_BOXPLOT = join("output", "{tax_level}-{method}-log2-{microbe}-reads-{celltype}-boxplot.png")

CONTAM_READ_PLOT = join("output", "{tax_level}-{method}_contaminant_microbe_reads_{celltype}.png")
NON_CONTAM_READ_PLOT = join("output", "{tax_level}-{method}_non_contaminant_microbe_reads_{celltype}.png")
MICROBE_READ_BOXPLOT = join("output", "{tax_level}-{method}-{microbe}_microbe_reads-{celltype}-boxplot.png")
MICROBE_READ_SWARMPLOT = join("output", "{patient}-{kingdom}-{tax_level}-{method}-{microbe}_microbe_reads-{celltype}-swarmplot.png")
MICROBE_CPM_PLOT = join("output", "{tax_level}-{method}-{microbe}_microbe_cpm-{celltype}.png")
MICROBE_CPM_SWARMPLOT = join("output", "{patient}-{kingdom}-{tax_level}-{method}-{microbe}-{celltype}-cpm-swarmplot.png")
NUM_HUMAN_READS_MICROBE_READS_PLOT = join("output", "{tax_level}-{method}_total_reads_human_microbes-{celltype}.png")
NUM_HUMAN_GENES_MICROBE_OTU_PLOT = join("output", "microbes_{tax_level}-{method}_human_genes-{celltype}.png")
DECONV_NORM_PLOT = join("output", "deconv-normalization-{patient}-{celltype}-{microbe}-{tax_level}-{kingdom}-{method}-plot.png")
ALL_DECONV_NORM_PLOT = join("output", "deconv-normalization-{celltype}-{microbe}-{tax_level}-{method}-plot.png")

HUMAN_SPIKE_RATIO_SWARMPLOT = join("output", "{patient}-{celltype}-human-ERCC-ratio-swarmplot.png")
MICROBE_HUMAN_SPIKE_RATIO_SWARMPLOT = join("output", "{patient}-{kingdom}-{tax_level}-{method}-{microbe}-{celltype}-human-ERCC-normalization-swarmplot.png")
# Visualization Rules
rule plot_deconv_normalization_all_patients:
    wildcard_constraints:
        celltype="patient"
    input:
        MICROBE_READ_TABLE,
        SAMPLE_METADATA,
    output:
        ALL_DECONV_NORM_PLOT
    script:
        "../src/plot_deconv_normalized_counts.R"

rule plot_deconv_normalization:
    input:
        PATIENT_MICROBE_READ_TABLE,
        PATIENT_SAMPLE_METADATA,
    output:
        DECONV_NORM_PLOT
    script:
        "../src/plot_deconv_normalized_counts.R"

rule plot_microbe_reads_by_human_reads:
    input:
        STAR_READCOUNT_TABLE,
        MICROBE_READ_TABLE,
        SAMPLE_METADATA
    output:
        NUM_HUMAN_READS_MICROBE_READS_PLOT,
        NUM_HUMAN_GENES_MICROBE_OTU_PLOT
    script:
        "../src/plot_human_reads_by_microbe_reads.py"

rule boxplot_microbe_reads_across_celltypes:
    input:
        MICROBE_READ_TABLE,
        SAMPLE_METADATA,
    output:
        LOG2_MICROBE_READ_BOXPLOT
    script:
        "../src/boxplot_log2_microbe_reads.py"


rule swarmplot_microbe_normalized_by_host_spikein_ratio_across_celltypes:
    input:
        STAR_READCOUNT_TABLE,
        MICROBE_READ_TABLE,
        SAMPLE_METADATA,
    output:
        MICROBE_HUMAN_SPIKE_RATIO_SWARMPLOT
    script:
        "../src/swarmplot_microbe_using_host_spikein_ratio.py"

rule swarmplot_host_spikein_ratio_across_celltypes:
    input:
        STAR_READCOUNT_TABLE,
        "output/genus_PathSeq_Bacteria_metadata.tsv",
    output:
        HUMAN_SPIKE_RATIO_SWARMPLOT
    script:
        "../src/swarmplot_host_spikein_ratio.py"

rule swarmplot_microbe_reads_across_celltypes:
    input:
        PATIENT_MICROBE_READ_TABLE,
        PATIENT_SAMPLE_METADATA,
    output:
        MICROBE_READ_SWARMPLOT
    script:
        "../src/swarmplot_microbe_readcounts.py"

rule swarmplot_microbe_cpms_across_celltypes:
    input:
        PATIENT_MICROBE_READ_TABLE,
        PATIENT_SAMPLE_METADATA,
    output:
        MICROBE_CPM_SWARMPLOT
    script:
        "../src/swarmplot_microbe_cpm.py"

rule boxplot_log2_totalreads_across_celltypes:
    input:
        MICROBE_READ_TABLE,
        SAMPLE_METADATA
    output:
        LOG2_TOTAL_READS_BOXPLOT
    script:
        "../src/boxplot_log2_total_reads.py"

rule boxplot_totalreads_across_celltypes:
    input:
        MICROBE_READ_TABLE,
        SAMPLE_METADATA
    output:
        TOTAL_READS_BOXPLOT
    script:
        "../src/boxplot_total_reads.py"

# rule plot_microbe_cpm_across_celltypes:
#     input:
#         MICROBE_CPM_TABLE,
#         SAMPLE_METADATA,
#     output:
#         MICROBE_CPM_PLOT
#     script:
#         "../src/plot_microbe_readcounts.py"
