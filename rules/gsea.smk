

PATIENT_HOST_GENE_MICROBE_ABUNDANCE_CORRELATION = join(
    "output", "{patient}", "host-gene-microbe-abundance-correlation",
    "spearman-correlation-{kingdom}-{tax_level}-{microbe}-{celltype}-{celltype_of_interest}-{celltype_comparison}-{block}-{method}.tsv"
    )
GSEA_CORR_TABLE = join(
    "output", "{patient}", "GSEA-microbe-abundance-correlation",
    "spearman-correlation-{kingdom}-{tax_level}-{microbe}-{celltype}-{celltype_of_interest}-{celltype_comparison}-{block}-{method}.tsv"
    )


rule run_GSEA_GO_terms_correlation:
    conda:
        "../envs/fgsea-env.yaml"
    input:
        "../data/c5.go.v7.2.symbols.gmt",
        PATIENT_HOST_GENE_MICROBE_ABUNDANCE_CORRELATION
    output:
        GSEA_CORR_TABLE
    script:
        "../src/run_fgsea_corr.R"


rule calculate_host_gene_microbe_abundance_correlation:
    wildcard_constraints:
        norm="spike"
    input:
        PATIENT_HUMAN_READCOUNT,
        PATIENT_SPIKEIN_READCOUNT,
        PATIENT_MICROBE_READ_TABLE,
        PATIENT_SAMPLE_METADATA,
        PATHSEQ_TAXID_MAP
    output:
        PATIENT_HOST_GENE_MICROBE_ABUNDANCE_CORRELATION
    script:
        "../src/compute_microbe_abundance_host_gene_expression_correlation.R"
