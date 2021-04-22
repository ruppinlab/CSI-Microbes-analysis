from os.path import join

# directories
STARsolo_PE_dir = join("raw", "star", "{patient}-{sample}-{plate}", "_STARpe", "Solo.out", "Gene", "raw")
STARsolo_SE_dir = join("raw", "star", "{patient}-{sample}-{plate}", "_STARse", "Solo.out", "Gene", "raw")

# STARsolo Input Files
STARsolo_PE_features = join(STARsolo_PE_dir, "features.tsv")
STARsolo_PE_barcodes = join(STARsolo_PE_dir, "barcodes.tsv")
STARsolo_PE_matrix = join(STARsolo_PE_dir, "matrix.mtx")
STARsolo_PE_features_gz = join(STARsolo_PE_dir, "features.tsv.gz")
STARsolo_PE_barcodes_gz = join(STARsolo_PE_dir, "barcodes.tsv.gz")
STARsolo_PE_matrix_gz = join(STARsolo_PE_dir, "matrix.mtx.gz")
STARsolo_SE_features = join(STARsolo_SE_dir, "features.tsv")
STARsolo_SE_barcodes = join(STARsolo_SE_dir, "barcodes.tsv")
STARsolo_SE_matrix = join(STARsolo_SE_dir, "matrix.mtx")
STARsolo_SE_features_gz = join(STARsolo_SE_dir, "features.tsv.gz")
STARsolo_SE_barcodes_gz = join(STARsolo_SE_dir, "barcodes.tsv.gz")
STARsolo_SE_matrix_gz = join(STARsolo_SE_dir, "matrix.mtx.gz")

# STARsolo Output Files
PLATE_STARsolo_READCOUNT = join("output", "star", "{patient}-{sample}-{plate}", "matrix.tsv")
PATIENT_STARsolo_READCOUNT = join("output", "star", "{patient}", "matrix.tsv")
SAMPLE_STARsolo_READCOUNT = join("output", "star", "{patient}-{sample}", "matrix.tsv")


def get_patient_GE_files(wildcards):
    p = units.loc[units.patient == wildcards.patient][["patient", "sample", "plate"]].drop_duplicates()
    return expand(PLATE_STARsolo_READCOUNT, zip, patient=p["patient"], sample=p["sample"], plate=p["plate"])

def get_sample_GE_files(wildcards):
    p = units.loc[units["sample"] == wildcards.sample][["patient", "sample", "plate"]].drop_duplicates()
    return expand(PLATE_STARsolo_READCOUNT, zip, patient=p["patient"], sample=p["sample"], plate=p["plate"])




rule combine_gene_expression_one_sample:
    input:
        get_sample_GE_files
    output:
        SAMPLE_STARsolo_READCOUNT
    script:
        "../src/concatenate_STARsolo_plate.py"

rule combine_gene_expression_one_patient:
    input:
        get_patient_GE_files
    output:
        PATIENT_STARsolo_READCOUNT
    script:
        "../src/concatenate_STARsolo_plate.py"

rule convert_STARsolo_output_to_csv:
    params:
        STARsolo_PE_dir,
        STARsolo_SE_dir
    input:
        STARsolo_PE_features_gz,
        STARsolo_PE_barcodes_gz,
        STARsolo_PE_matrix_gz,
        STARsolo_SE_features_gz,
        STARsolo_SE_barcodes_gz,
        STARsolo_SE_matrix_gz,
    output:
        PLATE_STARsolo_READCOUNT
    script:
        "../src/convert_STARsolo_output_to_tsv.R"

rule gzip_STARsolo_file:
    input:
        STARsolo_PE_features,
        STARsolo_PE_barcodes,
        STARsolo_PE_matrix,
        STARsolo_SE_features,
        STARsolo_SE_barcodes,
        STARsolo_SE_matrix
    output:
        STARsolo_PE_features_gz,
        STARsolo_PE_barcodes_gz,
        STARsolo_PE_matrix_gz,
        STARsolo_SE_features_gz,
        STARsolo_SE_barcodes_gz,
        STARsolo_SE_matrix_gz,
    shell:
        "gzip {input[0]} && "
        "gzip {input[1]} && "
        "gzip {input[2]} && "
        "gzip {input[3]} && "
        "gzip {input[4]} && "
        "gzip {input[5]}"
