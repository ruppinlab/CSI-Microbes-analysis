rule calculate_markers:
    wildcard_constraints:
        norm="deconv",
        kingdom="Bacteria"
    input:
        MICROBE_READ_TABLE,
        SAMPLE_METADATA
    output:
        TTEST_MARKERS,
        WILCOX_MARKERS
    script:
        "../src/run_scran_marker_analysis.R"

rule calculate_patient_markers:
    wildcard_constraints:
        norm="deconv",
        kingdom="Bacteria"
    input:
        PATIENT_MICROBE_READ_TABLE,
        PATIENT_SAMPLE_METADATA
    output:
        PATIENT_TTEST_MARKERS,
        PATIENT_WILCOX_MARKERS
    script:
        "../src/run_scran_marker_analysis.R"

rule calculate_sample_markers:
    wildcard_constraints:
        norm="deconv",
        kingdom="Bacteria"
    input:
        SAMPLE_MICROBE_READ_TABLE,
        SAMPLE_SAMPLE_METADATA
    output:
        SAMPLE_TTEST_MARKERS,
        SAMPLE_WILCOX_MARKERS
    script:
        "../src/run_scran_marker_analysis.R"

rule calculate_plate_markers:
    wildcard_constraints:
        norm="deconv",
        kingdom="Bacteria"
    input:
        PLATE_MICROBE_READ_TABLE,
        PLATE_SAMPLE_METADATA
    output:
        PLATE_TTEST_MARKERS,
        PLATE_WILCOX_MARKERS
    script:
        "../src/run_scran_marker_analysis.R"

rule calculate_patient_binomial_markers:
    input:
        PATIENT_MICROBE_READ_TABLE,
        PATIENT_SAMPLE_METADATA
    output:
        PATIENT_BINOM_MARKERS,
    script:
        "../src/run_scran_binomial_marker_analysis.R"
