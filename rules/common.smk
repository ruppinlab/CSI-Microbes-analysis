wildcard_constraints:
    tax_level="root|superkingdom|kingdom|phylum|class|order|family|genus|species|strain|no_rank",
    kingdom="Bacteria|Archaea|Viruses|Fungi|Eukaryota|root|All",
    method="PathSeq"


# patient-specific intermediate files
PATIENT_MICROBE_READ_TABLE = join("output", "{patient}", "{tax_level}_{method}_{kingdom}_reads.tsv")
PATIENT_SAMPLE_METADATA = join("output", "{patient}", "{method}_metadata.tsv")

SAMPLE_MICROBE_READ_TABLE = join("output", "{patient}", "{sample}", "{tax_level}_{method}_{kingdom}_reads.tsv")
#SAMPLE_SAMPLE_METADATA = join("output", "{patient}", "{sample}", "{method}_metadata.tsv")

PLATE_MICROBE_READ_TABLE = join("output", "{patient}", "{sample}", "{plate}", "{tax_level}_{method}_{kingdom}_reads.tsv")
#PLATE_SAMPLE_METADATA = join("output", "{patient}", "{sample}", "{plate}", "{method}_metadata.tsv")


### rules for splitting microbial read count files by plate or sample ###

rule split_sample_read_table_by_plate:
    input:
        SAMPLE_MICROBE_READ_TABLE,
        PATIENT_SAMPLE_METADATA
    output:
        PLATE_MICROBE_READ_TABLE,
    script:
        "../src/split_read_matrices_by_plate.py"

rule split_patient_read_table_by_sample:
    input:
        PATIENT_MICROBE_READ_TABLE,
        PATIENT_SAMPLE_METADATA
    output:
        SAMPLE_MICROBE_READ_TABLE,
    script:
        "../src/split_read_matrices_by_sample.py"
