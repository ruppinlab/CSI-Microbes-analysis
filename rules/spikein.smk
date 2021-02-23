

PATIENT_SPIKE_NORM_PLOT = join("output", "plots", "{patient}", "spike-normalization-{celltype}-{microbe}-{tax_level}-{kingdom}-{method}-plot.png")
SAMPLE_SPIKE_NORM_PLOT = join("output", "plots", "{patient}", "{sample}", "spike-normalization-{celltype}-{microbe}-{tax_level}-{kingdom}-{method}-plot.png")
PLATE_SPIKE_NORM_PLOT = join("output", "plots", "{patient}", "{sample}", "{plate}", "spike-normalization-{celltype}-{microbe}-{tax_level}-{kingdom}-{method}-plot.png")

PATIENT_SPIKEIN_READCOUNT = join("output", "star", "{patient}", "spike_ins.tsv")
PATIENT_HUMAN_READCOUNT = join("output", "star", "{patient}", "human_genes.tsv")
SAMPLE_SPIKEIN_READCOUNT = join("output", "star", "{patient}-{sample}", "spike_ins.tsv")
SAMPLE_HUMAN_READCOUNT = join("output", "star", "{patient}-{sample}", "human_genes.tsv")
PLATE_SPIKEIN_READCOUNT = join("output", "star", "{patient}-{sample}-{plate}", "spike_ins.tsv")
PLATE_HUMAN_READCOUNT = join("output", "star", "{patient}-{sample}-{plate}", "human_genes.tsv")

TTEST_MARKERS = join("output", "t-test-{celltype}-{celltype_of_interest}-{celltype_comparison}={tax_level}-{method}-{norm}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")
WILCOX_MARKERS = join("output", "wilcox-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{norm}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")

PATIENT_TTEST_MARKERS = join("output", "{patient}", "t-test-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{norm}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")
PATIENT_WILCOX_MARKERS = join("output", "{patient}", "wilcox-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{norm}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")

SAMPLE_TTEST_MARKERS = join("output", "{patient}", "{sample}", "t-test-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{norm}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")
SAMPLE_WILCOX_MARKERS = join("output", "{patient}", "{sample}", "wilcox-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{norm}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")

PLATE_TTEST_MARKERS = join("output", "{patient}", "{sample}", "{plate}", "t-test-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{norm}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")
PLATE_WILCOX_MARKERS = join("output", "{patient}", "{sample}", "{plate}", "wilcox-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{norm}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")

PATIENT_WILCOX_MARKERS_TAXA_NAME = join("output", "{patient}", "wilcox-taxa-names-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{norm}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")

rule split_patient_readcount:
    params:
        spike=SPIKE_PREFIX
    input:
        PATIENT_STARsolo_READCOUNT,
    output:
        PATIENT_SPIKEIN_READCOUNT,
        PATIENT_HUMAN_READCOUNT
    script:
        "../src/split_readcount.py"

rule split_sample_readcount:
    params:
        spike=SPIKE_PREFIX
    input:
        SAMPLE_STARsolo_READCOUNT,
    output:
        SAMPLE_SPIKEIN_READCOUNT,
        SAMPLE_HUMAN_READCOUNT
    script:
        "../src/split_readcount.py"

rule split_plate_readcount:
    params:
        spike=SPIKE_PREFIX
    input:
        PLATE_STARsolo_READCOUNT,
    output:
        PLATE_SPIKEIN_READCOUNT,
        PLATE_HUMAN_READCOUNT
    script:
        "../src/split_readcount.py"

rule calculate_patient_markers_using_spikeins:
    wildcard_constraints:
        norm="spike"
    params:
        spike_functions="../src/spike_normalization_functions.R"
    input:
        PATIENT_MICROBE_READ_TABLE,
        PATIENT_SPIKEIN_READCOUNT,
        PATIENT_SAMPLE_METADATA
    output:
        PATIENT_TTEST_MARKERS,
        PATIENT_WILCOX_MARKERS,
    script:
        "../src/run_scran_spikein_marker_analysis.R"

rule calculate_sample_markers_using_spikeins:
    wildcard_constraints:
        norm="spike"
    params:
        spike_functions="../src/spike_normalization_functions.R"
    input:
        SAMPLE_MICROBE_READ_TABLE,
        SAMPLE_SPIKEIN_READCOUNT,
        SAMPLE_SAMPLE_METADATA
    output:
        SAMPLE_TTEST_MARKERS,
        SAMPLE_WILCOX_MARKERS,
    script:
        "../src/run_scran_spikein_marker_analysis.R"

# rules to plot using spike-in normalization
rule plot_spike_normalization_patient:
    params:
        spike_functions="../src/spike_normalization_functions.R"
    input:
        PATIENT_MICROBE_READ_TABLE,
        PATIENT_SPIKEIN_READCOUNT,
        PATIENT_SAMPLE_METADATA,
        PATHSEQ_TAXID_MAP
    output:
        PATIENT_SPIKE_NORM_PLOT
    script:
        "../src/plot_spikein_normalized_counts.R"

rule plot_spike_normalization_sample:
    params:
        spike_functions="../src/spike_normalization_functions.R"
    input:
        SAMPLE_MICROBE_READ_TABLE,
        SAMPLE_SPIKEIN_READCOUNT,
        SAMPLE_SAMPLE_METADATA,
        PATHSEQ_TAXID_MAP
    output:
        SAMPLE_SPIKE_NORM_PLOT
    script:
        "../src/plot_spikein_normalized_counts.R"

rule plot_spike_normalization_plate:
    params:
        spike_functions="../src/spike_normalization_functions.R"
    input:
        PLATE_MICROBE_READ_TABLE,
        PLATE_SPIKEIN_READCOUNT,
        PLATE_SAMPLE_METADATA,
        PATHSEQ_TAXID_MAP
    output:
        PLATE_SPIKE_NORM_PLOT
    script:
        "../src/plot_spikein_normalized_counts.R"

rule convert_taxid_to_name:
    input:
        PATIENT_WILCOX_MARKERS,
        PATHSEQ_TAXID_MAP
    output:
        PATIENT_WILCOX_MARKERS_TAXA_NAME
    script:
        "../src/convert_taxid_to_taxa_name.R"
