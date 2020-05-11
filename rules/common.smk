wildcard_constraints:
    tax_level="root|superkingdom|phylum|class|order|family|genus|species"


SAMPLES_FILE = join("data", "samples.tsv")
#ruleorder: convert_patient_to_CPM > convert_patient_batch_to_CPM
# contamination files
CONTAMINANTS_FILE = join("output", "contaminants.tsv")
POORE2020_SUPP_TABLE = join("..", "data", "Poore2020_SuppTables.xlsx")

# read count file
STAR_READCOUNT_TABLE = join("output", "star-readcounts.tsv")

# intermediate files
MICROBE_READ_TABLE = join("output", "{tax_level}_{method}_microbe_reads.tsv")
MICROBE_CPM_TABLE = join("output", "{tax_level}_{method}_microbe_cpm.tsv")
SAMPLE_METADATA = join("output", "{tax_level}_{method}_sample_metadata.tsv")


# patient-specific intermediate files
PATIENT_MICROBE_READ_TABLE = join("output", "{patient}_{tax_level}_{method}_microbe_reads.tsv")
PATIENT_SAMPLE_METADATA = join("output", "{patient}_{tax_level}_{method}_metadata.tsv")
#PATIENT_BATCH_MICROBE_READ_TABLE = join("output", "{patient}_{tax_level}_microbe_{batch}_reads.tsv")
#PATIENT_BATCH_SAMPLE_METADATA = join("output", "{patient}_{tax_level}_{batch}_metadata.tsv")
PATIENT_MICROBE_CPM_TABLE = join("output", "{patient}_{tax_level}_{method}_microbe_cpm.tsv")
#PATIENT_BATCH_MICROBE_CPM_TABLE = join("output", "{patient}_{tax_level}_{batch}_cpm.tsv")
# output files

# TOP_MICROBES_FILE = join("output", "{tax_level}_microbes_ranked.tsv")

# rule convert_Kraken_to_counts:
#     wildcard_constraints:
#         tax_level="genus",
#         method="Kraken"
#     input:
#         SAMPLES_FILE,
#         #paired=expand(KRAKEN_PAIRED_BIOM, zip, patient=samples["patient"], sample=samples["sample"]),
#         #unpaired=expand(KRAKEN_UNPAIRED_BIOM, zip, patient=samples["patient"], sample=samples["sample"]),
#     output:
#         MICROBE_READ_TABLE,
#         SAMPLE_METADATA
#     script:
#         "../src/convert_Kraken_to_counts.py"
#
# # rules for contaminants
# rule clean_contaminants_file:
#     input:
#         POORE2020_SUPP_TABLE
#     output:
#         CONTAMINANTS_FILE
#     script:
#         "../src/clean_contaminants.py"

# rules for processing data
rule filter_read_table_by_patient:
    input:
        MICROBE_READ_TABLE,
        SAMPLE_METADATA
    output:
        PATIENT_MICROBE_READ_TABLE,
        PATIENT_SAMPLE_METADATA
    script:
        "../src/split_read_matrices_by_patient.py"

# rule filter_patient_read_table_by_batch:
#     input:
#         PATIENT_MICROBE_READ_TABLE,
#         PATIENT_SAMPLE_METADATA
#     output:
#         PATIENT_BATCH_MICROBE_READ_TABLE,
#         PATIENT_BATCH_SAMPLE_METADATA
#     script:
#         "../src/split_read_matrices_by_batch.py"

# rule convert_patient_batch_to_CPM:
#     input:
#         PATIENT_BATCH_MICROBE_READ_TABLE
#     output:
#         PATIENT_BATCH_MICROBE_CPM_TABLE
#     script:
#         "../src/convert_to_CPM.py"

# rule convert_patient_to_CPM:
#     input:
#         PATIENT_MICROBE_READ_TABLE
#     output:
#         PATIENT_MICROBE_CPM_TABLE
#     script:
#         "../src/convert_to_CPM.py"
#
# rule convert_to_CPM:
#     input:
#         MICROBE_READ_TABLE
#     output:
#         MICROBE_CPM_TABLE
#     script:
#         "../src/convert_to_CPM.py"
