# CSI-Microbes-analysis

## Setting up the environment

```
conda env create -f envs/CSI-Microbes.yaml
conda activate CSI-Microbes
```

## Generating Aulicino2018 Results

First, you will need the download the STAR and PathSeq files for Aulicino2018 and then you can run snakemake.

```
cd Aulicino2018

rsync -avc --include='pathseq.txt' --include='*/' --exclude='*' helix:/data/Robinson-SB/scRNA-seq-microbe-identification/Aulicino2018/output/ data/
rsync -avc --include='barcodes.tsv' --include='*/' --exclude='*' helix:/data/Robinson-SB/scRNA-seq-microbe-identification/Aulicino2018/output/ data/
rsync -avc --include='features.tsv' --include='*/' --exclude='*' helix:/data/Robinson-SB/scRNA-seq-microbe-identification/Aulicino2018/output/ data/
rsync -avc --include='matrix.mtx' --include='*/' --exclude='*' helix:/data/Robinson-SB/scRNA-seq-microbe-identification/Aulicino2018/output/ data/

snakemake --cores 2 all
```

## Generating Maynard2020 Results

First, you will need the download the STAR and PathSeq files for Maynard2020 (which can take a while) and then you can run snakemake.

```
cd Maynard2020

rsync -avc --include='pathseq.txt' --include='*/' --exclude='*' helix:/data/Robinson-SB/scRNA-seq-microbe-identification/Maynard2020/output/ data/

rsync -avc --include='barcodes.tsv' --include='*/' --exclude='*' helix:/data/Robinson-SB/scRNA-seq-microbe-identification/Maynard2020/output/ data/
rsync -avc --include='features.tsv' --include='*/' --exclude='*' helix:/data/Robinson-SB/scRNA-seq-microbe-identification/Maynard2020/output/ data/
rsync -avc --include='matrix.mtx' --include='*/' --exclude='*' helix:/data/Robinson-SB/scRNA-seq-microbe-identification/Maynard2020/output/ data/

snakemake --cores 1 all
```
