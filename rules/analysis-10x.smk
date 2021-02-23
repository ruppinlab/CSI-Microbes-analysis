### rules for analyzing microbial reads from 10x ###

PATIENT_MICROBE_READ_TABLE = join("output", "{patient}", "{tax_level}_{method}_{kingdom}_reads.tsv")
PATIENT_SAMPLE_METADATA = join("output", "{patient}", "{tax_level}_{method}_{kingdom}_metadata.tsv")
SAMPLE_MICROBE_READ_TABLE = join("output", "{patient}", "{sample}_{tax_level}_{method}_{kingdom}_reads.tsv")
SAMPLE_SAMPLE_METADATA = join("output", "{patient}", "{sample}_{tax_level}_{method}_{kingdom}_metadata.tsv")

BINOM_MARKERS = join("output", "binomial-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")
PATIENT_BINOM_MARKERS = join("output", "{patient}", "binomial-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")
SAMPLE_BINOM_MARKERS = join("output", "{patient}", "{sample}", "binomial-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")
PATIENT_BINOM_MARKERS_TAXA_NAME = join("output", "{patient}", "binomial-taxa-names-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")

rule plot_sample_binomial_result:
    input:
        SAMPLE_MICROBE_READ_TABLE,
        SAMPLE_SAMPLE_METADATA
    output:
        SAMPLE_BINOM_PLOT
    script:
        "../src/plot_binomial_result.R"

### rules to calculate the binomials at the sample or patient level ###
rule calculate_sample_binomial_markers:
    input:
        SAMPLE_MICROBE_READ_TABLE,
        SAMPLE_SAMPLE_METADATA
    output:
        SAMPLE_BINOM_MARKERS,
    script:
        "../src/run_scran_binomial_marker_analysis.R"

rule calculate_patient_binomial_markers:
    input:
        PATIENT_MICROBE_READ_TABLE,
        PATIENT_SAMPLE_METADATA
    output:
        PATIENT_BINOM_MARKERS,
    script:
        "../src/run_scran_binomial_marker_analysis.R"

rule convert_binomial_markers_taxid_to_name:
    input:
        PATIENT_BINOM_MARKERS,
        PATHSEQ_TAXID_MAP
    output:
        PATIENT_BINOM_MARKERS_TAXA_NAME
    script:
        "../src/convert_taxid_to_taxa_name.R"

### convert PathSeq output files to a matrix of microbial counts ###

rule convert_PathSeq_to_readcounts:
    wildcard_constraints:
        method="PathSeq"
    params:
        join("data", "PathSeq", "{}-{}-{}", "pathseq.txt")
    input:
        "data/units.tsv"
    output:
        PATIENT_MICROBE_READ_TABLE,
        PATIENT_SAMPLE_METADATA
    script:
        "../src/convert_10x-PathSeq_output_to_read_counts.py"
