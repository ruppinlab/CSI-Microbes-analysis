
# patient-specific intermediate files
PATIENT_MICROBE_READ_TABLE = join("output", "{patient}", "{tax_level}_{method}_{kingdom}_reads.tsv")
PATIENT_SAMPLE_METADATA = join("output", "{patient}", "{tax_level}_{method}_{kingdom}_metadata.tsv")

PATIENT_PATHSEQ_EDGELIST_FILE = join("output", "{patient}", "edgelist_{kingdom}_{method}.tsv")
PATIENT_PATHSEQ_TAXID_MAP = join("output", "{patient}", "tax_id_map_{kingdom}_{method}.tsv")


rule extract_PathSeq_edgelist:
    params:
        join("raw", "PathSeq", "{}-{}-{}-{}", "pathseq.txt")
    input:
        "data/units.tsv"
    output:
        PATIENT_PATHSEQ_EDGELIST_FILE
    script:
        "../src/extract_SS2_PathSeq_edgelist.py"

rule extract_name_tax_id_mapping:
    params:
        join("raw", "PathSeq", "{}-{}-{}-{}", "pathseq.txt")
    input:
        "data/units.tsv"
    output:
        PATIENT_PATHSEQ_TAXID_MAP
    script:
        "../src/extract_SS2_PathSeq_name_tax_id.py"

rule convert_PathSeq_to_read_counts:
    wildcard_constraints:
        method="PathSeq"
    params:
        join("raw", "PathSeq", "{}-{}-{}-{}", "pathseq.txt")
    input:
        "data/units.tsv",
    output:
        PATIENT_MICROBE_READ_TABLE,
        PATIENT_SAMPLE_METADATA,
    script:
        "../src/convert_SS2-PathSeq_output_to_read_counts.py"
