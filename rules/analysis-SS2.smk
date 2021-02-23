
# patient-specific intermediate files
PATIENT_MICROBE_READ_TABLE = join("output", "{patient}", "{tax_level}_{method}_{kingdom}_reads.tsv")
PATIENT_SAMPLE_METADATA = join("output", "{patient}", "{tax_level}_{method}_{kingdom}_metadata.tsv")


rule convert_PathSeq_to_read_counts:
    wildcard_constraints:
        method="PathSeq"
    params:
        join("data", "PathSeq", "{}-{}-{}-{}", "pathseq.txt")
    input:
        "data/units.tsv",
    output:
        PATIENT_MICROBE_READ_TABLE,
        PATIENT_SAMPLE_METADATA,
    script:
        "../src/convert_SS2-PathSeq_output_to_read_counts.py"
