# CSI-Microbes-analysis

## Setting up the environment

```
conda env create -f envs/CSI-Microbes.yaml
conda activate CSI-Microbes-env
```

## Ben-Moshe2019 Figures

First, you will need the PathSeq and SRPRISM output files for Ben-Moshe2019. To download these data files, you will need to be connected to the NIH network and have access to the `/data/Robinson-SB` directory and follow the below instructions

```
cd Ben-Moshe2019
rsync -avc --include='pathseq.txt' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Ben-Moshe2019/output/ data/
rsync -avc --include='CB-UMI-count-SL1344.tsv' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Ben-Moshe2019/output/ data/
```

Next, you can generate the figures using the below command

```
snakemake --cores <number of CPUs> plot_figures
```

## Aulicino2018 Figures

First, you will need the PathSeq, SRPRISM and STAR output files for Aulicino2018. To download these data files, you will need to be connected to the NIH network and have access to the `/data/Robinson-SB` directory and follow the below instructions

```
cd Aulicino2018
rsync -avc --include='barcodes.tsv' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Aulicino2018/output/ data/
rsync -avc --include='features.tsv' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Aulicino2018/output/ data/
rsync -avc --include='matrix.mtx' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Aulicino2018/output/ data/
rsync -avc --include='pathseq.txt' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Aulicino2018/output/ data/
rsync -avc --include='*-paired-count.gff' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Aulicino2018/output/ data/
```

Next, you can generate the figures using the below command

```
snakemake --cores <number of CPUs> plot_figures
```

## Lee2020 Figures

First, you will need the PathSeq output files for Lee2020. To download these data files, you will need to be connected to the NIH network and have access to the `/data/Robinson-SB` directory and follow the below instructions

```
cd Lee2020
rsync -avc --include='pathseq.txt' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Lee2020/output/ data/
```

Next, you can generate the figures using the below command

```
snakemake --cores <number of CPUs> plot_figures
```

## Maynard2020 Figures

First, you will need the PathSeq output files for Maynard2020. To download these data files, you will need to be connected to the NIH network and have access to the `/data/Robinson-SB` directory and follow the below instructions

```
cd Maynard2020
rsync -avc --include='pathseq.txt' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Maynard2020/output/ data/
rsync -avc --include='barcodes.tsv' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Maynard2020/output/ data/
rsync -avc --include='features.tsv' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Maynard2020/output/ data/
rsync -avc --include='matrix.mtx' --include='*/' --exclude='*' helix:/data/Robinson-SB/CSI-Microbes-identification/Maynard2020/output/ data/
```

Next, you can generate the figures using the below command

```
snakemake --cores <number of CPUs> plot_figures
```
