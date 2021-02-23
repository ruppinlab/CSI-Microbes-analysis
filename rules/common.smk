wildcard_constraints:
    tax_level="root|superkingdom|kingdom|phylum|class|order|family|genus|species|strain|no_rank",
    kingdom="Bacteria|Archaea|Viruses|Fungi|Eukaryota|root|All",
    method="PathSeq"


# patient-specific intermediate files
PATIENT_MICROBE_READ_TABLE = join("output", "{patient}", "{tax_level}_{method}_{kingdom}_reads.tsv")
PATIENT_SAMPLE_METADATA = join("output", "{patient}", "{tax_level}_{method}_{kingdom}_metadata.tsv")

SAMPLE_MICROBE_READ_TABLE = join("output", "{patient}", "{sample}_{tax_level}_{method}_{kingdom}_reads.tsv")
SAMPLE_SAMPLE_METADATA = join("output", "{patient}", "{sample}_{tax_level}_{method}_{kingdom}_metadata.tsv")

PLATE_MICROBE_READ_TABLE = join("output", "{patient}", "{sample}_{plate}_{tax_level}_{method}_{kingdom}_reads.tsv")
PLATE_SAMPLE_METADATA = join("output", "{patient}", "{sample}_{plate}_{tax_level}_{method}_{kingdom}_metadata.tsv")

PATHSEQ_EDGELIST_FILE = join("output", "{patient}", "edgelist_{kingdom}_{method}.tsv")
PATHSEQ_TAXID_MAP = join("output", "{patient}", "tax_id_map_{kingdom}_{method}.tsv")



rule extract_PathSeq_edgelist:
    params:
        join("data", "PathSeq", "{}-{}-{}-{}", "pathseq.txt")
    input:
        "data/units.tsv"
    output:
        PATHSEQ_EDGELIST_FILE
    script:
        "../src/extract_PathSeq_edgelist.py"

rule extract_name_tax_id_mapping:
    params:
        join("data", "PathSeq", "{}-{}-{}-{}", "pathseq.txt")
    input:
        "data/units.tsv"
    output:
        PATHSEQ_TAXID_MAP
    script:
        "../src/extract_PathSeq_name_tax_id.py"

### rules for splitting microbial read count files by plate or sample ###

rule split_sample_metadata_by_plate:
    input:
        SAMPLE_MICROBE_READ_TABLE,
        SAMPLE_SAMPLE_METADATA
    output:
        PLATE_MICROBE_READ_TABLE,
        PLATE_SAMPLE_METADATA
    script:
        "../src/split_read_matrices_by_plate.py"

rule split_patient_metadata_by_sample:
    input:
        PATIENT_MICROBE_READ_TABLE,
        PATIENT_SAMPLE_METADATA
    output:
        SAMPLE_MICROBE_READ_TABLE,
        SAMPLE_SAMPLE_METADATA
    script:
        "../src/split_read_matrices_by_sample.py"
