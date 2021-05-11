

rule convert_to_CPM_human_reads:
    input:
        MICROBE_READ_TABLE,
        STAR_READCOUNT_TABLE
    output:
        MICROBE_CPM_TABLE
    script:
        "../src/convert_to_CPM_human_reads.py"

rule convert_reads_to_CPM_per_patient:
    input:
        PATIENT_MICROBE_READ_TABLE,
        STAR_READCOUNT_TABLE
    output:
        PATIENT_MICROBE_CPM_TABLE
    script:
        "../src/convert_to_CPM.py"

rule plot_microbe_cpm_across_celltypes:
    input:
        MICROBE_CPM_TABLE,
        MICROBE_SAMPLE_METADATA,
    output:
        MICROBE_CPM_PLOT
    script:
        "../src/plot_microbe_readcounts.py"

rule plot_contaminant_reads_across_celltypes:
    input:
        MICROBE_READ_TABLE,
        MICROBE_SAMPLE_METADATA,
        POORE2020_SUPP_TABLE
    output:
        CONTAM_READ_PLOT,
        NON_CONTAM_READ_PLOT
    script:
        "../src/plot_contaminant_readcounts.py"

rule plot_contaminant_cpm_across_celltypes:
    input:
        MICROBE_CPM_TABLE,
        MICROBE_SAMPLE_METADATA,
        POORE2020_SUPP_TABLE
    output:
        CONTAM_CPM_PLOT,
        NON_CONTAM_CPM_PLOT
    script:
        "../src/plot_contaminant_readcounts.py"

rule remove_likely_contaminants:
    wildcard_constraints:
        tax_level="genus"
    input:
        MICROBE_READ_TABLE,
        POORE2020_SUPP_TABLE
    output:
        NO_CONTAM_MICROBE_READ_TABLE
    script:
        "../src/remove_contaminants.py"

rule identify_most_common_microbes:
    input:
        MICROBE_READ_TABLE,
        CONTAMINANTS_FILE
    output:
        TOP_MICROBES_FILE
    script:
        "../src/identify_top_microbes.py"

rule compute_wilcoxon_rank_sum:
    input:
        MICROBE_READ_TABLE,
        MICROBE_SAMPLE_METADATA,
        CONTAMINANTS_FILE
    output:
        RANKSUM_RESULT_FILE
    script:
        "../src/compute_wilcoxon_rank_sum_test.py"

rule run_edgeR_DE_human_seqdepth:
    input:
        PATIENT_MICROBE_READ_TABLE,
        PATIENT_MICROBE_SAMPLE_METADATA,
        STAR_READCOUNT_TABLE
    output:
        EDGER_HUMAN_SEQDEPTH_RESULTS
    script:
        "../src/run_edgeR_human_seqdepth.R"
