
wildcard_constraints:
    seed="\d+"

READ_SUBSAMPLED_PATIENT_SAMPLE_METADATA = join("read-subsampling-output", "{patient}_{celltype}_{celltype_of_interest}_{celltype_comparison}_{tax_level}_{method}_{kingdom}_{microbe}_{nreads}_{seed}_metadata.tsv")
READ_SUBSAMPLED_PATIENT_MICROBE_READ_TABLE = join("read-subsampling-output", "{patient}_{celltype}_{celltype_of_interest}_{celltype_comparison}_{tax_level}_{method}_{kingdom}_{microbe}_{nreads}_{seed}_reads.tsv")
READ_SUBSAMPLED_PATIENT_HUMAN_MICROBE_READ_TABLE = join("read-subsampling-output", "{patient}_{celltype}_{celltype_of_interest}_{celltype_comparison}_{tax_level}_{method}_{kingdom}_{microbe}_{nreads}_{seed}_and_human_reads.tsv")

READ_SUBSAMPLED_TTEST_MARKERS = join("read-subsampling-output", "t-test-{patient}_{celltype}_{celltype_of_interest}_{celltype_comparison}_{tax_level}_{method}_{kingdom}_{microbe}_{nreads}_{seed}_{norm}_{pvaltype}_{lfc}_{block}.tsv")
READ_SUBSAMPLED_WILCOX_MARKERS = join("read-subsampling-output", "wilcox-{patient}_{celltype}_{celltype_of_interest}_{celltype_comparison}_{tax_level}_{method}_{kingdom}_{microbe}_{nreads}_{seed}_{norm}_{pvaltype}_{lfc}_{block}.tsv")


READ_SUBSAMPLED_WILCOX_SIG_RES = join("read-subsampling-output", "sigpercentage-wilcox-{patient}_{celltype}_{celltype_of_interest}_{celltype_comparison}_{tax_level}_{method}_{microbe}_{kingdom}_{nsamples}_{nreads}_{norm}_{pvaltype}_{lfc}_{block}.tsv")
READ_SAMPLE_SIZE_PLOT = join("output", "sigpercentage-wilcox-{patient}_{celltype}_{celltype_of_interest}_{celltype_comparison}_{tax_level}_{method}_{microbe}_{kingdom}_{nsamples}_{minreads}_{maxreads}_{norm}_{pvaltype}_{lfc}_{block}.pdf")


def get_sig_files(wildcards):
    minreads = int(wildcards["minreads"])
    maxreads = int(wildcards["maxreads"])
    return expand(READ_SUBSAMPLED_WILCOX_SIG_RES, patient="{patient}", celltype="{celltype}",
           celltype_of_interest="{celltype_of_interest}", tax_level="{tax_level}", microbe="{microbe}",
           method="{method}", kingdom="{kingdom}", norm="{norm}", celltype_comparison="{celltype_comparison}",
           pvaltype="{pvaltype}", nsamples="{nsamples}", lfc="{lfc}",
           block="{block}", nreads=range(minreads, maxreads, 100))

rule plot_sample_size_analysis:
    input:
        get_sig_files
    output:
        READ_SAMPLE_SIZE_PLOT
    script:
        "../src/plot_sample_size_analysis.py"

def get_wilcox_input_files(wildcards):
    nsamples = int(wildcards["nsamples"])
    return expand(READ_SUBSAMPLED_WILCOX_MARKERS, patient="{patient}", celltype="{celltype}",
           celltype_of_interest="{celltype_of_interest}",
           celltype_comparison="{celltype_comparison}", tax_level="{tax_level}",
           method="{method}", kingdom="{kingdom}", nreads="{nreads}",
           norm="{norm}", pvaltype="{pvaltype}", lfc="{lfc}", block="{block}",
           microbe="{microbe}", seed=range(0,nsamples))

# {patient}_{celltype}_{celltype_of_interest}_{tax_level}_{method}_{kingdom}_{ncells}_{seed}{norm}_{pvaltype}
rule calculate_num_significant:
    input:
        get_wilcox_input_files
    output:
        READ_SUBSAMPLED_WILCOX_SIG_RES
    script:
        "../src/count_significant_results.py"

rule compute_deconv_significance:
    wildcard_constraints:
        norm="deconv"
    input:
        READ_SUBSAMPLED_PATIENT_MICROBE_READ_TABLE,
        READ_SUBSAMPLED_PATIENT_SAMPLE_METADATA,
    output:
        READ_SUBSAMPLED_TTEST_MARKERS,
        READ_SUBSAMPLED_WILCOX_MARKERS
    script:
        "../src/run_scran_marker_analysis.R"

rule compute_spike_significance:
    wildcard_constraints:
        norm="spike"
    params:
        spike=SPIKE_PREFIX
    input:
        READ_SUBSAMPLED_PATIENT_HUMAN_MICROBE_READ_TABLE,
        READ_SUBSAMPLED_PATIENT_SAMPLE_METADATA,
    output:
        READ_SUBSAMPLED_TTEST_MARKERS,
        READ_SUBSAMPLED_WILCOX_MARKERS
    script:
        "../src/run_scran_spikein_marker_analysis.R"

rule combine_subsampled_microbial_human_reads:
    input:
        STAR_READCOUNT_TABLE,
        READ_SUBSAMPLED_PATIENT_MICROBE_READ_TABLE
    output:
        READ_SUBSAMPLED_PATIENT_HUMAN_MICROBE_READ_TABLE
    script:
        "../src/combine_microbial_human_reads.py"

rule downsample_microbial_reads:
    input:
        PATIENT_SAMPLE_METADATA,
        PATIENT_MICROBE_READ_TABLE
    output:
        READ_SUBSAMPLED_PATIENT_SAMPLE_METADATA,
        READ_SUBSAMPLED_PATIENT_MICROBE_READ_TABLE
    script:
        "../src/downsample_microbial_reads.py"
