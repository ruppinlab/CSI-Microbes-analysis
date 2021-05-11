
PATIENT_MICROBE_PLATE_ENRICHMENT_FILE = join("output", "{batch}-enrichment-{tax_level}-{patient}.tsv")

rule identify_microbes_associated_with_plate_per_patient:
    input:
        PATIENT_MICROBE_READ_TABLE,
        PATIENT_SAMPLE_METADATA
    output:
        PATIENT_MICROBE_PLATE_ENRICHMENT_FILE
    script:
        "../src/identify_microbes_associated_with_plate.py"
