# Aulicino2018

## Setting up the environment

```
conda env create -f ../envs/CSI-Microbes.yaml
conda activate CSI-Microbes
```

## Generating Figure 1

First, you will need the PathSeq output files, which should go in `data/PathSeq`. For analyses using spike-in normalization, you will need to STAR readcount files, which should go in `data/STAR`.

```
snakemake --cores 1
```
