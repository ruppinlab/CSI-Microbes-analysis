### rules for analyzing microbial reads from 10x ###

PATIENT_MICROBE_READ_TABLE = join("output", "{patient}", "{tax_level}_{method}_{kingdom}_reads.tsv")
PATIENT_SAMPLE_METADATA = join("output", "{patient}", "{tax_level}_{method}_{kingdom}_metadata.tsv")
SAMPLE_MICROBE_READ_TABLE = join("output", "{patient}", "{sample}", "{tax_level}_{method}_{kingdom}_reads.tsv")
SAMPLE_SAMPLE_METADATA = join("output", "{patient}", "{sample}", "{tax_level}_{method}_{kingdom}_metadata.tsv")

# BINOM_MARKERS = join("output", "binomial-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")
# PATIENT_BINOM_MARKERS = join("output", "{patient}", "binomial-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")
# PATIENT_BINOM_VOLCANO = join("output", "plots", "{patient}", "volcano-binomial-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.png")
# SAMPLE_BINOM_MARKERS = join("output", "{patient}", "{sample}", "binomial-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")
# SAMPLE_BINOM_VOLCANO = join("output", "plots", "{patient}", "{sample}", "volcano-binomial-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.png")
# SAMPLE_BINOM_PLOT = join("output", "{patient}", "{sample}", "plots", "binomial-plot-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{kingdom}-{pvaltype}-{lfc}-{block}-{microbe_of_interest}.jpg")

SAMPLE_FISHER_EXACT = join("output", "{patient}", "{sample}", "fisher-exact-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{kingdom}-{minprop}.tsv")


PATIENT_PATHSEQ_EDGELIST_FILE = join("output", "{patient}", "edgelist_{kingdom}_{method}.tsv")
PATIENT_PATHSEQ_TAXID_MAP = join("output", "{patient}", "tax_id_map_{kingdom}_{method}.tsv")

rule extract_PathSeq_edgelist:
    params:
        join("raw", "PathSeq", "{}-{}-{}", "pathseq.txt")
    input:
        "data/units.tsv"
    output:
        PATIENT_PATHSEQ_EDGELIST_FILE
    script:
        "../src/extract_10x_PathSeq_edgelist.py"

rule extract_name_tax_id_mapping:
    params:
        join("raw", "PathSeq", "{}-{}-{}", "pathseq.txt")
    input:
        "data/units.tsv"
    output:
        PATIENT_PATHSEQ_TAXID_MAP
    script:
        "../src/extract_10x_PathSeq_name_tax_id.py"

rule calculate_Fishers_exact_sample:
    input:
        SAMPLE_MICROBE_READ_TABLE,
        SAMPLE_SAMPLE_METADATA,
        PATIENT_PATHSEQ_TAXID_MAP
    output:
        SAMPLE_FISHER_EXACT
    script:
        "../src/calculate_Fisher_exact_test.R"

# rule plot_sample_binomial_result:
#     input:
#         SAMPLE_MICROBE_READ_TABLE,
#         SAMPLE_SAMPLE_METADATA,
#         PATIENT_PATHSEQ_TAXID_MAP
#     output:
#         SAMPLE_BINOM_PLOT
#     script:
#         "../src/plot_binomial_result.R"

### rules to calculate the binomials at the sample or patient level ###
# rule calculate_sample_binomial_markers:
#     input:
#         SAMPLE_MICROBE_READ_TABLE,
#         SAMPLE_SAMPLE_METADATA,
#         PATIENT_PATHSEQ_TAXID_MAP
#     output:
#         SAMPLE_BINOM_MARKERS,
#         #SAMPLE_BINOM_VOLCANO
#     script:
#         "../src/run_scran_binomial_marker_analysis.R"
#
# rule calculate_patient_binomial_markers:
#     input:
#         PATIENT_MICROBE_READ_TABLE,
#         PATIENT_SAMPLE_METADATA,
#         PATIENT_PATHSEQ_TAXID_MAP
#     output:
#         PATIENT_BINOM_MARKERS,
#         #PATIENT_BINOM_VOLCANO
#     script:
#         "../src/run_scran_binomial_marker_analysis.R"

# rule convert_binomial_markers_taxid_to_name:
#     input:
#         PATIENT_BINOM_MARKERS,
#         PATHSEQ_TAXID_MAP
#     output:
#         PATIENT_BINOM_MARKERS_TAXA_NAME
#     script:
#         "../src/convert_taxid_to_taxa_name.R"

### convert PathSeq output files to a matrix of microbial counts ###

rule convert_PathSeq_to_readcounts:
    wildcard_constraints:
        method="PathSeq"
    params:
        join("raw", "PathSeq", "{}-{}-{}", "pathseq.txt")
    input:
        "data/units.tsv"
    output:
        PATIENT_MICROBE_READ_TABLE,
        PATIENT_SAMPLE_METADATA
    script:
        "../src/convert_10x-PathSeq_output_to_read_counts.py"
