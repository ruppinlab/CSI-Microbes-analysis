

PATIENT_RANKSUM_RESULT_FILE = join("output", "ranksum_{patient}_{tax_level}_{celltype}.tsv")
PATIENT_HYPERGEOM_RESULT_FILE = join("output", "hypergeometric_{patient}_{tax_level}_{celltype}.tsv")
PATIENT_BATCH_RANKSUM_RESULT_FILE = join("output", "ranksum_{patient}_{batch}_{tax_level}_{celltype}.tsv")
#RANKSUM_RESULT_FILE = join("output", "{tax_level}_ranksum_result-{celltype}.tsv")
#PATIENT_RANKSUM_CONTAMINANT_PLOT = join("output", "{tax_level}_ranksum_contaminants_{patient}-{celltype}_plot.png")

rule calculate_markers:
    wildcard_constraints:
        norm="deconv"
    input:
        PATIENT_MICROBE_READ_TABLE,
        PATIENT_SAMPLE_METADATA
    output:
        TTEST_MARKERS,
        WILCOX_MARKERS
    script:
        "../src/run_scran_marker_analysis.R"

# rules for running edgeR
# rule run_edgeR_DE:
#     input:
#         PATIENT_MICROBE_READ_TABLE,
#         PATIENT_SAMPLE_METADATA,
#     output:
#         EDGER_RESULTS
#     script:
#         "../src/run_edgeR_analysis.R"
#
# # rules for statistical analyses
# rule compute_wilcoxon_rank_sum_patient_species:
#     wildcard_constraints:
#         tax_level="species"
#     input:
#         PATIENT_MICROBE_CPM_TABLE,
#         PATIENT_SAMPLE_METADATA
#     output:
#         PATIENT_RANKSUM_RESULT_FILE
#     script:
#         "../src/compute_wilcoxon_rank_sum_test.py"
#
# rule compute_wilcoxon_rank_sum_patient_genus:
#     wildcard_constraints:
#         tax_level="genus"
#     input:
#         PATIENT_MICROBE_CPM_TABLE,
#         PATIENT_SAMPLE_METADATA,
#         CONTAMINANTS_FILE
#     output:
#         PATIENT_RANKSUM_RESULT_FILE
#     script:
#         "../src/compute_wilcoxon_rank_sum_test_contamination.py"

# rule compute_wilcoxon_rank_sum_patient_batch_genus:
#     wildcard_constraints:
#         tax_level="genus"
#     input:
#         PATIENT_BATCH_MICROBE_CPM_TABLE,
#         PATIENT_BATCH_SAMPLE_METADATA,
#         CONTAMINANTS_FILE
#     output:
#         PATIENT_BATCH_RANKSUM_RESULT_FILE
#     script:
#         "../src/compute_wilcoxon_rank_sum_test_contamination.py"

# hypergeometric enrichment test rules

# rule compute_hypergeom_patient_genus:
#     wildcard_constraints:
#         tax_level="genus"
#     input:
#         PATIENT_MICROBE_CPM_TABLE,
#         PATIENT_SAMPLE_METADATA,
#         CONTAMINANTS_FILE
#     output:
#         PATIENT_HYPERGEOM_RESULT_FILE
#     script:
#         "../src/compute_hypergeometric_test_contamination.py"
#
# rule compute_hypergeom_patient_species:
#     wildcard_constraints:
#         tax_level="species"
#     input:
#         PATIENT_MICROBE_CPM_TABLE,
#         PATIENT_SAMPLE_METADATA,
#         CONTAMINANTS_FILE
#     output:
#         PATIENT_HYPERGEOM_RESULT_FILE
#     script:
#         "../src/compute_hypergeometric_test.py"
# rule compute_contaminant_rank_sum_test:
#     input:
#         PATIENT_RANKSUM_RESULT_FILE
#     output:
#         PATIENT_RANKSUM_CONTAMINANT_PLOT
#     script:
#         "../src/plot_rank_sum_contaminants.py"
