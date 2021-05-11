include: "Snakefile"

FIGURE_S2 = join("output", "plots", "figure_S2.pdf")

SRPRISM_READ_COUNT_FILE = join("raw", "SRPRISM", "{patient}", "{sample}", "CB-UMI-count-{genome}.tsv")


rule plot_figure_S2:
    input:
        SAMPLE_SAMPLE_METADATA.format(patient="Pt0", sample="GSM3454529", tax_level="genus", method="PathSeq", kingdom="Bacteria"),
        SRPRISM_READ_COUNT_FILE.format(patient="Pt0", sample="GSM3454529", genome="SL1344")
    output:
        FIGURE_S2
    script:
        "src/plot_figure_S2.R"
