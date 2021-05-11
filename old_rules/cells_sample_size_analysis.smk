

SUBSAMPLED_PATIENT_SAMPLE_METADATA = join("output", "{patient}_{celltype}_{celltype_of_interest}_{celltype_comparison}_{tax_level}_{method}_{kingdom}_{ncells}_{seed}_metadata.tsv")
SUBSAMPLED_PATIENT_MICROBE_READ_TABLE = join("output", "{patient}_{celltype}_{celltype_of_interest}_{celltype_comparison}_{tax_level}_{method}_{kingdom}_{ncells}_{seed}_reads.tsv")

SUBSAMPLED_TTEST_MARKERS = join("output", "t-test-{patient}_{celltype}_{celltype_of_interest}_{celltype_comparison}_{tax_level}_{method}_{kingdom}_{ncells}_{seed}_{norm}_{pvaltype}.tsv")
SUBSAMPLED_WILCOX_MARKERS = join("output", "wilcox-{patient}_{celltype}_{celltype_of_interest}_{celltype_comparison}_{tax_level}_{method}_{kingdom}_{ncells}_{seed}_{norm}_{pvaltype}.tsv")


SUBSAMPLED_WILCOX_SIG_RES = join("output", "sigpercentage-wilcox-{patient}_{celltype}_{celltype_of_interest}_{celltype_comparison}_{tax_level}_{method}_{microbe}_{kingdom}_{nsamples}_{ncells}_{norm}_{pvaltype}.tsv")
SAMPLE_SIZE_PLOT = join("output", "sigpercentage-wilcox-{patient}_{celltype}_{celltype_of_interest}_{celltype_comparison}_{tax_level}_{method}_{microbe}_{kingdom}_{nsamples}_{mincells}_{maxcells}_{norm}_{pvaltype}.pdf")


def get_sig_files(wildcards):
    mincells = int(wildcards["mincells"])
    maxcells = int(wildcards["maxcells"])
    return expand(SUBSAMPLED_WILCOX_SIG_RES, patient="{patient}", celltype="{celltype}",
           celltype_of_interest="{celltype_of_interest}", tax_level="{tax_level}", microbe="{microbe}",
           method="{method}", kingdom="{kingdom}", norm="{norm}", celltype_comparison="{celltype_comparison}",
           pvaltype="{pvaltype}", nsamples="{nsamples}", ncells=range(mincells, maxcells, 5))

rule plot_sample_size_analysis:
    input:
        get_sig_files
    output:
        SAMPLE_SIZE_PLOT
    script:
        "../src/plot_sample_size_analysis.py"

def get_wilcox_input_files(wildcards):
    nsamples = int(wildcards["nsamples"])
    return expand(SUBSAMPLED_WILCOX_MARKERS, patient="{patient}", celltype="{celltype}",
           celltype_of_interest="{celltype_of_interest}",
           celltype_comparison="{celltype_comparison}", tax_level="{tax_level}",
           method="{method}", kingdom="{kingdom}", ncells="{ncells}",
           norm="{norm}", pvaltype="all", seed=range(0,nsamples))

# {patient}_{celltype}_{celltype_of_interest}_{tax_level}_{method}_{kingdom}_{ncells}_{seed}{norm}_{pvaltype}
rule calculate_num_significant:
    input:
        get_wilcox_input_files
    output:
        SUBSAMPLED_WILCOX_SIG_RES
    script:
        "../src/count_significant_results.py"

rule compute_significance:
    input:
        SUBSAMPLED_PATIENT_MICROBE_READ_TABLE,
        SUBSAMPLED_PATIENT_SAMPLE_METADATA,
    output:
        SUBSAMPLED_TTEST_MARKERS,
        SUBSAMPLED_WILCOX_MARKERS
    script:
        "../src/run_scran_marker_analysis.R"

rule subsample_cells:
    input:
        PATIENT_SAMPLE_METADATA,
        PATIENT_MICROBE_READ_TABLE
    output:
        SUBSAMPLED_PATIENT_SAMPLE_METADATA,
        SUBSAMPLED_PATIENT_MICROBE_READ_TABLE
    script:
        "../src/subsample_cells.py"
