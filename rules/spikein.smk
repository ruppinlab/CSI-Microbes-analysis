
ALL_MICROBE_SPIKEIN_CORR_FILE = join("output", "all-samples-{method}-{tax_level}-microbe-spikein-correlation.tsv")
PATIENT_MICROBE_SPIKEIN_CORR_FILE = join("output", "{patient}-{method}-{tax_level}-{celltype}-microbe-spikein-correlation.tsv")
ALL_SAMPLE_MICROBE_SPIKEIN_SCATTER = join("output", "all-samples-{method}-{tax_level}-{microbe}-spikein-scatterplot.png")


rule calculate_markers_using_spikeins:
    wildcard_constraints:
        norm="spike"
    params:
        spike=SPIKE_PREFIX
    input:
        PATIENT_MICROBE_READ_TABLE,
        PATIENT_SAMPLE_METADATA,
        STAR_READCOUNT_TABLE
    output:
        TTEST_MARKERS,
        WILCOX_MARKERS
    script:
        "../src/run_scran_spikein_marker_analysis.R"

rule plot_spike_normalization:
    params:
        spike=SPIKE_PREFIX
    input:
        PATIENT_MICROBE_READ_TABLE,
        PATIENT_SAMPLE_METADATA,
        STAR_READCOUNT_TABLE
    output:
        SPIKE_NORM_PLOT
    script:
        "../src/plot_spike_normalized_counts.R"

rule plot_normalization_approaches:
    params:
        spike=SPIKE_PREFIX
    input:
        PATIENT_MICROBE_READ_TABLE,
        PATIENT_SAMPLE_METADATA,
        STAR_READCOUNT_TABLE
    output:
        NORM_COMP_PLOT
    script:
        "../src/plot_normalization_approaches.R"

rule scatterplot_microbe_spikein_all_samples:
    params:
        spike=SPIKE_PREFIX
    input:
        STAR_READCOUNT_TABLE,
        MICROBE_READ_TABLE,
        SAMPLE_METADATA
    output:
        ALL_SAMPLE_MICROBE_SPIKEIN_SCATTER
    script:
        "../src/scatterplot_spikein_microbe.py"

rule calculate_microbe_spikein_correlation_patient:
    params:
        spike=SPIKE_PREFIX
    input:
        STAR_READCOUNT_TABLE,
        PATIENT_MICROBE_READ_TABLE,
        PATIENT_SAMPLE_METADATA
    output:
        PATIENT_MICROBE_SPIKEIN_CORR_FILE
    script:
        "../src/correlate_microbe_spikein_patient.py"

rule calculate_microbe_spikein_correlation_across_all_samples:
    params:
        spike=SPIKE_PREFIX
    input:
        STAR_READCOUNT_TABLE,
        MICROBE_READ_TABLE
    output:
        ALL_MICROBE_SPIKEIN_CORR_FILE
    script:
        "../src/identify_microbes_correlated_spike_ins.py"

# rule run_edgeR_spike_in_normalization:
#     params:
#         spike=SPIKE_PREFIX
#     input:
#         PATIENT_MICROBE_READ_TABLE,
#         PATIENT_SAMPLE_METADATA,
#         STAR_READCOUNT_TABLE
#     output:
#         EDGER_SPIKE_RESULTS
#     script:
#         "../src/run_edgeR_human_spikein.R"
