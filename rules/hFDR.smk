
#ruleorder: calculate_sample_binomial_markers > calculate_patient_binomial_markers

PATIENT_TTEST_MARKERS = join("output", "{patient}", "t-test-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{norm}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")
PATIENT_WILCOX_MARKERS = join("output", "{patient}", "wilcox-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{norm}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")
PATIENT_BINOM_MARKERS = join("output", "{patient}", "binomial-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")

SAMPLE_TTEST_MARKERS = join("output", "{patient}", "{sample}", "t-test-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{norm}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")
SAMPLE_WILCOX_MARKERS = join("output", "{patient}", "{sample}", "wilcox-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{norm}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")
SAMPLE_BINOM_MARKERS = join("output", "{patient}", "{sample}", "binomial-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")
SAMPLE_FISHER_MARKERS = join("output", "{patient}", "{sample}", "fisher-exact-{celltype}-{celltype_of_interest}-{celltype_comparison}-{tax_level}-{method}-{kingdom}-{minprop}.tsv")


# Wilcox files
CLASS_PATIENT_WILCOX_MARKERS = PATIENT_WILCOX_MARKERS.format(
    tax_level="class", kingdom="{kingdom}", method="{method}", patient="{patient}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", norm="{norm}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )
ORDER_PATIENT_WILCOX_MARKERS = PATIENT_WILCOX_MARKERS.format(
    tax_level="order", kingdom="{kingdom}", method="{method}", patient="{patient}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", norm="{norm}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )
FAMILY_PATIENT_WILCOX_MARKERS = PATIENT_WILCOX_MARKERS.format(
    tax_level="family", kingdom="{kingdom}", method="{method}", patient="{patient}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", norm="{norm}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )
GENUS_PATIENT_WILCOX_MARKERS = PATIENT_WILCOX_MARKERS.format(
    tax_level="genus", kingdom="{kingdom}", method="{method}", patient="{patient}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", norm="{norm}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )
SPECIES_PATIENT_WILCOX_MARKERS = PATIENT_WILCOX_MARKERS.format(
    tax_level="species", kingdom="{kingdom}", method="{method}", patient="{patient}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", norm="{norm}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )

CLASS_SAMPLE_WILCOX_MARKERS = SAMPLE_WILCOX_MARKERS.format(
    tax_level="class", kingdom="{kingdom}", method="{method}", patient="{patient}", sample="{sample}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", norm="{norm}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )
ORDER_SAMPLE_WILCOX_MARKERS = SAMPLE_WILCOX_MARKERS.format(
    tax_level="order", kingdom="{kingdom}", method="{method}", patient="{patient}", sample="{sample}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", norm="{norm}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )
FAMILY_SAMPLE_WILCOX_MARKERS = SAMPLE_WILCOX_MARKERS.format(
    tax_level="family", kingdom="{kingdom}", method="{method}", patient="{patient}", sample="{sample}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", norm="{norm}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )
GENUS_SAMPLE_WILCOX_MARKERS = SAMPLE_WILCOX_MARKERS.format(
    tax_level="genus", kingdom="{kingdom}", method="{method}", patient="{patient}", sample="{sample}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", norm="{norm}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )
SPECIES_SAMPLE_WILCOX_MARKERS = SAMPLE_WILCOX_MARKERS.format(
    tax_level="species", kingdom="{kingdom}", method="{method}", patient="{patient}", sample="{sample}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", norm="{norm}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )

# T-test files
CLASS_PATIENT_TTEST_MARKERS = PATIENT_TTEST_MARKERS.format(
    tax_level="class", kingdom="{kingdom}", method="{method}", patient="{patient}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", norm="{norm}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )
ORDER_PATIENT_TTEST_MARKERS = PATIENT_TTEST_MARKERS.format(
    tax_level="order", kingdom="{kingdom}", method="{method}", patient="{patient}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", norm="{norm}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )
FAMILY_PATIENT_TTEST_MARKERS = PATIENT_TTEST_MARKERS.format(
    tax_level="family", kingdom="{kingdom}", method="{method}", patient="{patient}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", norm="{norm}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )
GENUS_PATIENT_TTEST_MARKERS = PATIENT_TTEST_MARKERS.format(
    tax_level="genus", kingdom="{kingdom}", method="{method}", patient="{patient}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", norm="{norm}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )
SPECIES_PATIENT_TTEST_MARKERS = PATIENT_TTEST_MARKERS.format(
    tax_level="species", kingdom="{kingdom}", method="{method}", patient="{patient}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", norm="{norm}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )

CLASS_SAMPLE_TTEST_MARKERS = SAMPLE_TTEST_MARKERS.format(
    tax_level="class", kingdom="{kingdom}", method="{method}", patient="{patient}", sample="{sample}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", norm="{norm}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )
ORDER_SAMPLE_TTEST_MARKERS = SAMPLE_TTEST_MARKERS.format(
    tax_level="order", kingdom="{kingdom}", method="{method}", patient="{patient}", sample="{sample}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", norm="{norm}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )
FAMILY_SAMPLE_TTEST_MARKERS = SAMPLE_TTEST_MARKERS.format(
    tax_level="family", kingdom="{kingdom}", method="{method}", patient="{patient}", sample="{sample}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", norm="{norm}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )
GENUS_SAMPLE_TTEST_MARKERS = SAMPLE_TTEST_MARKERS.format(
    tax_level="genus", kingdom="{kingdom}", method="{method}", patient="{patient}", sample="{sample}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", norm="{norm}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )
SPECIES_SAMPLE_TTEST_MARKERS = SAMPLE_TTEST_MARKERS.format(
    tax_level="species", kingdom="{kingdom}", method="{method}", patient="{patient}", sample="{sample}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", norm="{norm}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )

# Binomial test files
CLASS_PATIENT_BINOM_MARKERS = PATIENT_BINOM_MARKERS.format(
    tax_level="class", kingdom="{kingdom}", method="{method}", patient="{patient}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )
ORDER_PATIENT_BINOM_MARKERS = PATIENT_BINOM_MARKERS.format(
    tax_level="order", kingdom="{kingdom}", method="{method}", patient="{patient}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )
FAMILY_PATIENT_BINOM_MARKERS = PATIENT_BINOM_MARKERS.format(
    tax_level="family", kingdom="{kingdom}", method="{method}", patient="{patient}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )
GENUS_PATIENT_BINOM_MARKERS = PATIENT_BINOM_MARKERS.format(
    tax_level="genus", kingdom="{kingdom}", method="{method}", patient="{patient}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )
SPECIES_PATIENT_BINOM_MARKERS = PATIENT_BINOM_MARKERS.format(
    tax_level="species", kingdom="{kingdom}", method="{method}", patient="{patient}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )

CLASS_SAMPLE_BINOM_MARKERS = SAMPLE_BINOM_MARKERS.format(
    tax_level="class", kingdom="{kingdom}", method="{method}", patient="{patient}", sample="{sample}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )
ORDER_SAMPLE_BINOM_MARKERS = SAMPLE_BINOM_MARKERS.format(
    tax_level="order", kingdom="{kingdom}", method="{method}", patient="{patient}", sample="{sample}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )
FAMILY_SAMPLE_BINOM_MARKERS = SAMPLE_BINOM_MARKERS.format(
    tax_level="family", kingdom="{kingdom}", method="{method}", patient="{patient}", sample="{sample}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )
GENUS_SAMPLE_BINOM_MARKERS = SAMPLE_BINOM_MARKERS.format(
    tax_level="genus", kingdom="{kingdom}", method="{method}", patient="{patient}", sample="{sample}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )
SPECIES_SAMPLE_BINOM_MARKERS = SAMPLE_BINOM_MARKERS.format(
    tax_level="species", kingdom="{kingdom}", method="{method}", patient="{patient}", sample="{sample}", celltype="{celltype}", lfc="{lfc}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", pvaltype="{pvaltype}", block="{block}", direction="{direction}", minprop="{minprop}"
    )

# Fisher's exact files
CLASS_SAMPLE_FISHER_MARKERS = SAMPLE_FISHER_MARKERS.format(
    tax_level="class", kingdom="{kingdom}", method="{method}", patient="{patient}", sample="{sample}", celltype="{celltype}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", minprop="{minprop}"
    )
ORDER_SAMPLE_FISHER_MARKERS = SAMPLE_FISHER_MARKERS.format(
    tax_level="order", kingdom="{kingdom}", method="{method}", patient="{patient}", sample="{sample}", celltype="{celltype}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", minprop="{minprop}"
    )
FAMILY_SAMPLE_FISHER_MARKERS = SAMPLE_FISHER_MARKERS.format(
    tax_level="family", kingdom="{kingdom}", method="{method}", patient="{patient}", sample="{sample}", celltype="{celltype}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", minprop="{minprop}"
    )
GENUS_SAMPLE_FISHER_MARKERS = SAMPLE_FISHER_MARKERS.format(
    tax_level="genus", kingdom="{kingdom}", method="{method}", patient="{patient}", sample="{sample}", celltype="{celltype}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", minprop="{minprop}"
    )
SPECIES_SAMPLE_FISHER_MARKERS = SAMPLE_FISHER_MARKERS.format(
    tax_level="species", kingdom="{kingdom}", method="{method}", patient="{patient}", sample="{sample}", celltype="{celltype}", celltype_of_interest="{celltype_of_interest}", celltype_comparison="{celltype_comparison}", minprop="{minprop}"
    )

PATIENT_hFDR_WILCOX_MARKERS = join("output", "results", "{patient}",
                               "hFDR-wilcox-{celltype}-{celltype_of_interest}-{celltype_comparison}-{method}-{norm}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")

SAMPLE_hFDR_WILCOX_MARKERS = join("output", "results", "{patient}", "{sample}",
                               "hFDR-wilcox-{celltype}-{celltype_of_interest}-{celltype_comparison}-{method}-{norm}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")

PATIENT_hFDR_TTEST_MARKERS = join("output", "results", "{patient}",
                               "hFDR-ttest-{celltype}-{celltype_of_interest}-{celltype_comparison}-{method}-{norm}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")

SAMPLE_hFDR_TTEST_MARKERS = join("output", "results", "{patient}", "{sample}",
                               "hFDR-ttest-{celltype}-{celltype_of_interest}-{celltype_comparison}-{method}-{norm}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")

PATIENT_hFDR_BINOM_MARKERS = join("output", "results", "{patient}",
                               "hFDR-binom-{celltype}-{celltype_of_interest}-{celltype_comparison}-{method}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")

SAMPLE_hFDR_BINOM_MARKERS = join("output", "results", "{patient}", "{sample}",
                               "hFDR-binom-{celltype}-{celltype_of_interest}-{celltype_comparison}-{method}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")

SAMPLE_hFDR_FISHER_MARKERS = join("output", "results", "{patient}", "{sample}",
                               "hFDR-fisher-exact-{celltype}-{celltype_of_interest}-{celltype_comparison}-{method}-{kingdom}-{minprop}.tsv")


rule adjust_hFRD_wilcox_patient:
    params:
        metric="AUC"
    input:
        edgelist = PATIENT_PATHSEQ_EDGELIST_FILE,
        tax_id_map = PATIENT_PATHSEQ_TAXID_MAP,
        class_markers = CLASS_PATIENT_WILCOX_MARKERS,
        order_markers = ORDER_PATIENT_WILCOX_MARKERS,
        family_markers = FAMILY_PATIENT_WILCOX_MARKERS,
        genus_markers = GENUS_PATIENT_WILCOX_MARKERS,
        species_markers = SPECIES_PATIENT_WILCOX_MARKERS,
    output:
        PATIENT_hFDR_WILCOX_MARKERS
    script:
        "../src/run_hFDR.py"

rule adjust_hFRD_wilcox_sample:
    params:
        metric="AUC"
    input:
        edgelist = PATIENT_PATHSEQ_EDGELIST_FILE,
        tax_id_map = PATIENT_PATHSEQ_TAXID_MAP,
        class_markers = CLASS_SAMPLE_WILCOX_MARKERS,
        order_markers = ORDER_SAMPLE_WILCOX_MARKERS,
        family_markers = FAMILY_SAMPLE_WILCOX_MARKERS,
        genus_markers = GENUS_SAMPLE_WILCOX_MARKERS,
        species_markers = SPECIES_SAMPLE_WILCOX_MARKERS,
    output:
        SAMPLE_hFDR_WILCOX_MARKERS
    script:
        "../src/run_hFDR.py"

rule adjust_hFRD_ttest_patient:
    params:
        metric="logFC"
    input:
        edgelist = PATIENT_PATHSEQ_EDGELIST_FILE,
        tax_id_map = PATIENT_PATHSEQ_TAXID_MAP,
        class_markers = CLASS_PATIENT_TTEST_MARKERS,
        order_markers = ORDER_PATIENT_TTEST_MARKERS,
        family_markers = FAMILY_PATIENT_TTEST_MARKERS,
        genus_markers = GENUS_PATIENT_TTEST_MARKERS,
        species_markers = SPECIES_PATIENT_TTEST_MARKERS,
    output:
        PATIENT_hFDR_TTEST_MARKERS
    script:
        "../src/run_hFDR.py"

rule adjust_hFRD_ttest_sample:
    params:
        metric="logFC"
    input:
        edgelist = PATIENT_PATHSEQ_EDGELIST_FILE,
        tax_id_map = PATIENT_PATHSEQ_TAXID_MAP,
        class_markers = CLASS_SAMPLE_TTEST_MARKERS,
        order_markers = ORDER_SAMPLE_TTEST_MARKERS,
        family_markers = FAMILY_SAMPLE_TTEST_MARKERS,
        genus_markers = GENUS_SAMPLE_TTEST_MARKERS,
        species_markers = SPECIES_SAMPLE_TTEST_MARKERS,
    output:
        SAMPLE_hFDR_TTEST_MARKERS
    script:
        "../src/run_hFDR.py"

rule adjust_hFRD_binom_patient:
    params:
        metric="logFC"
    input:
        edgelist = PATIENT_PATHSEQ_EDGELIST_FILE,
        tax_id_map = PATIENT_PATHSEQ_TAXID_MAP,
        class_markers = CLASS_PATIENT_BINOM_MARKERS,
        order_markers = ORDER_PATIENT_BINOM_MARKERS,
        family_markers = FAMILY_PATIENT_BINOM_MARKERS,
        genus_markers = GENUS_PATIENT_BINOM_MARKERS,
        species_markers = SPECIES_PATIENT_BINOM_MARKERS,
    output:
        PATIENT_hFDR_BINOM_MARKERS
    script:
        "../src/run_hFDR.py"

rule adjust_hFRD_binom_sample:
    params:
        metric="logFC"
    input:
        edgelist = PATIENT_PATHSEQ_EDGELIST_FILE,
        tax_id_map = PATIENT_PATHSEQ_TAXID_MAP,
        class_markers = CLASS_SAMPLE_BINOM_MARKERS,
        order_markers = ORDER_SAMPLE_BINOM_MARKERS,
        family_markers = FAMILY_SAMPLE_BINOM_MARKERS,
        genus_markers = GENUS_SAMPLE_BINOM_MARKERS,
        species_markers = SPECIES_SAMPLE_BINOM_MARKERS,
    output:
        SAMPLE_hFDR_BINOM_MARKERS
    script:
        "../src/run_hFDR.py"

rule adjust_hFRD_fisher_sample:
    params:
        metric="odds.ratio"
    input:
        edgelist = PATIENT_PATHSEQ_EDGELIST_FILE,
        tax_id_map = PATIENT_PATHSEQ_TAXID_MAP,
        class_markers = CLASS_SAMPLE_FISHER_MARKERS,
        order_markers = ORDER_SAMPLE_FISHER_MARKERS,
        family_markers = FAMILY_SAMPLE_FISHER_MARKERS,
        genus_markers = GENUS_SAMPLE_FISHER_MARKERS,
        species_markers = SPECIES_SAMPLE_FISHER_MARKERS,
    output:
        SAMPLE_hFDR_FISHER_MARKERS
    script:
        "../src/run_hFDR.py"
