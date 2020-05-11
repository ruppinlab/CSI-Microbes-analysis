# CSI-Microbes-analysis

## Setting up the environment

```
conda env create -f envs/CSI-Microbes.yaml
conda activate CSI-Microbes
```

## Generating Figure 1

First, you will need the PathSeq output files for Aulicino2018, which should go in `Aulicino2018/data/PathSeq`. For analyses using spike-in normalization, you will need to STAR readcount files, which should go in `Aulicino2018/data/STAR`.

```
cd Aulicino2018
snakemake --cores 1 plot_fig1
```
