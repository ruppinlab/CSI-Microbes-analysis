# decontam files
DECONTAM_READ_TABLE = join("output", "decontam", "{tax_level}_microbe_read_table.tsv")
DECONTAM_METADATA = join("output", "decontam", "{tax_level}_metadata.tsv")
DECONTAM_NEGATIVE_CONTAMINANT_FILE = join("output", "decontam", "negative_control_contaminant_scores.tsv")

# run decontam
rule run_decontam_negative_controls:
    conda:
        "../envs/decontam-env.yaml"
    input:
        DECONTAM_READ_TABLE,
        DECONTAM_METADATA
    output:
        DECONTAM_NEGATIVE_CONTAMINANT_FILE
    script:
        "../src/run_decontam_negative_control.R"

rule identify_negative_control_samples:
    input:
        STAR_READCOUNT_TABLE,
        MICROBE_SAMPLE_METADATA,
        MICROBE_READ_TABLE
    output:
        DECONTAM_READ_TABLE,
        DECONTAM_METADATA
    script:
        "../src/identify_negative_control_samples.py"
