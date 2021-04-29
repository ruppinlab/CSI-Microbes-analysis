# CSI-Microbes-analysis

This repository contains part of the workflows for reproducing the results from the bioarxiv paper [Identifying the Landscape of Intratumoral Microbes via a Single Cell Transcriptomic Analysis](https://www.biorxiv.org/content/10.1101/2020.05.14.096230v1) by Welles Robinson, Fiorella Schischlik, Michael Gertz, Alejandro Schaffer and Eytan Ruppin.  This repository contains the workflows to analyze microbial reads from 10x and Smart-seq2 scRNA-seq datasets to identify microbial taxa that are differentially abundant or differentially present. Prior to running this code, these microbial reads must be identified using the [CSI-Microbes-identification repository](https://github.com/ruppinlab/CSI-Microbes-identification). The code in this repository was written by Welles Robinson and alpha-tested by Alejandro Schaffer.

## Setting up the environment

This workflow has minimal computational constraints (the major computational steps are in the CSI-Microbes-identification pipeline). I run these steps locally on my Mac with 32G memory without a problem. For now, this workflow requires that you are on the NIH network (either physically present or connected via the VPN) and have access to the biowulf directory `/data/Robinson-SB/CSI-Microbes-identification` to download the necessary files. 

This workflow expects that conda has been installed. For instructions on how to install conda, see [conda install documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

Next, the GitHub repository needs to be cloned as below. The below instructions assume that you have an ssh key associated with your GitHub account. If you do not, you can generate a new ssh key and associate it with your GitHub username by following [these instructions](https://docs.github.com/en/github/authenticating-to-github/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent).

```
git clone git@github.com:ruppinlab/CSI-Microbes-analysis.git
```

Next, you need to create the conda environment (you need to perform this step only once unless you explicitly delete the conda environment).

```
cd CSI-Microbes-analysis
conda env create -f envs/CSI-Microbes-analysis.yaml
```

Finally, you need to activate the recently created conda environment (all of the commands assume that the conda environment `CSI-Microbes-env` is active).

```
conda activate CSI-Microbes-env
```

## Software Dependencies

TODO

## Reproducing results from Aulicino2018

To reproduce the results from Aulicino2018, you will need to download the PathSeq, STAR and SRPRISM files using the below command (within the Aulicino2018 directory)

```
rsync -avc --include='pathseq.txt' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Aulicino2018/output/ raw/
rsync -avc --include='barcodes.tsv' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Aulicino2018/output/ raw/
rsync -avc --include='features.tsv' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Aulicino2018/output/ raw/
rsync -avc --include='matrix.mtx' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Aulicino2018/output/ raw/
rsync -avc --include='*-paired-count.gff' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Aulicino2018/output/ raw/
```

### Reproducing Figure 2A

To generate figure 2A (`output/plots/figure_2A_1.pdf` and `output/plots/figure_2A_2.pdf`) use the below command (within the Aulicino2018 directory)

```
snakemake --cores <number of CPUs> plot_figure2A
```

### Reproducing Figure S1

To generate figure S1 (`output/plots/figure_S1A.pdf`, `output/plots/figure_S1B.pdf`, `output/plots/figure_S1C.pdf`) use the below command (within the Aulicino2018 directory)

```
snakemake --cores <number of CPUs> plot_figure_S1
```

## Reproducing results from Ben-Moshe2019

To reproduce the results from Ben-Moshe2019, you will need to download the PathSeq and SRPRISM files using the below command (within the Ben-Moshe2019 directory)

```
rsync -avc --include='pathseq.txt' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Ben-Moshe2019/output/ raw/
rsync -avc --include='CB-UMI-count-SL1344.tsv' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Ben-Moshe2019/output/ raw/
```

### Reproducing Figure 2B

Next, you can generate figure 2B (`output/plots/figure_2B_1.pdf` and `output/plots/figure_2B_2.pdf`) using the below command

```
snakemake --cores <number of CPUs> plot_figure_2B
```

### Reproducing Figure S2

Next, you can generate figure S2 (`output/plots/figure_S2.pdf`) using the below command

```
snakemake --cores <number of CPUs> plot_figure_S2
```


## Reproducing results from Lee2020

To reproduce figure 3, you will need to download the PathSeq files using the below command (within the Lee2020 directory)

```
rsync -avc --include='pathseq.txt' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Lee2020/output/ raw/
```


Next, you can generate figure 3 (`output/plots/figure_3A_1.pdf`, `output/plots/figure_3A_2.pdf`, `output/plots/figure_3B_1.pdf` and `output/plots/figure_3B_2.pdf`) using the below command

```
snakemake --cores <number of CPUs> plot_figure3
```

<!-- ## Reproducing results from Paulson2018 -->

## Reproducing results from Maynard2020

To reproduce figure 4 and 5, you will need to download the PathSeq and STAR files using the below commands (within the Maynard2020 directory).

```
rsync -avc --include='pathseq.txt' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Maynard2020/output/ raw/
rsync -avc --include='barcodes.tsv' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Maynard2020/output/ raw/
rsync -avc --include='features.tsv' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Maynard2020/output/ raw/
rsync -avc --include='matrix.mtx' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Maynard2020/output/ raw/
```

### Reproducing Figure 4

Next, you can generate figure 4 (`output/plots/figure_4A_1.pdf`, `output/plots/figure_4A_2.pdf`, `output/plots/figure_4A_3.pdf`, `output/plots/figure_4B_1.pdf`, `output/plots/figure_4B_2.pdf`, `output/plots/figure_4B_3.pdf`, `output/plots/figure_4C_1.pdf`, `output/plots/figure_4C_2.pdf`) using the below command

```
snakemake --cores <number of CPUs> plot_figure_4
```

### Reproducing Figure 5A

Next, you can plot figure 5A (`output/plots/figure_5A.pdf`) using the below command

```
snakemake --cores <number of CPUs> plot_figure_5a
```
