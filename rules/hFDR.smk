
#ruleorder: calculate_sample_binomial_markers > calculate_patient_binomial_markers


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


PATIENT_hFDR_WILCOX_MARKERS = join("output", "results", "{patient}",
                               "hFDR-wilcox-{celltype}-{celltype_of_interest}-{celltype_comparison}-{method}-{norm}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")

SAMPLE_hFDR_WILCOX_MARKERS = join("output", "results", "{patient}", "{sample}",
                               "hFDR-wilcox-{celltype}-{celltype_of_interest}-{celltype_comparison}-{method}-{norm}-{kingdom}-{pvaltype}-{lfc}-{block}-{direction}-{minprop}.tsv")



rule adjust_hFRD_patient:
    input:
        edgelist = PATHSEQ_EDGELIST_FILE,
        tax_id_map = PATHSEQ_TAXID_MAP,
        class_markers = CLASS_PATIENT_WILCOX_MARKERS,
        order_markers = ORDER_PATIENT_WILCOX_MARKERS,
        family_markers = FAMILY_PATIENT_WILCOX_MARKERS,
        genus_markers = GENUS_PATIENT_WILCOX_MARKERS,
        species_markers = SPECIES_PATIENT_WILCOX_MARKERS,
    output:
        PATIENT_hFDR_WILCOX_MARKERS
    script:
        "../src/run_hFDR.py"

rule adjust_hFRD_sample:
    input:
        edgelist = PATHSEQ_EDGELIST_FILE,
        tax_id_map = PATHSEQ_TAXID_MAP,
        class_markers = CLASS_SAMPLE_WILCOX_MARKERS,
        order_markers = ORDER_SAMPLE_WILCOX_MARKERS,
        family_markers = FAMILY_SAMPLE_WILCOX_MARKERS,
        genus_markers = GENUS_SAMPLE_WILCOX_MARKERS,
        species_markers = SPECIES_SAMPLE_WILCOX_MARKERS,
    output:
        SAMPLE_hFDR_WILCOX_MARKERS
    script:
        "../src/run_hFDR.py"
