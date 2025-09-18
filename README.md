# SYMPHONY Screen SYM101
This repository contains the SYMPHONY Screen SYM101 pipeline, which is a Snakemake pipeline for processing the SYM101 CRISPR screen. It includes steps for data processing, statistical analysis, and visualization.

# Running this pipeline
- Step1: install conda / mamba if necessary (mamba is much faster)
- Step2: Clone this repository to your machine using `git clone`
- Step3: Go to `configuration/config.yaml` and change the `TEMPDIR` to the location where you want to store the temporary files. This will be quite a lot of data, so make sure you have enough space.
- Step4: Install snakemake using `conda install bioconda::snakemake=8.14.0` or `mamba install bioconda::snakemake=8.14.0`
- Step5: run go to the snakemake folder and run `snakemake --sdm conda` to run the pipeline

# What this pipeline does
![Snakemake rulegraph of the SYM101 pipeline](snakemake/snake_rulegraph.png?raw=true "SYM101 pipeline")

- Process raw data from the sequence reads to count matrix
- Run statistical analysis tools (BAGEL, MAGeCK, gscreend) to produce result tables
- Figures are then plotted in notebooks (see `notebooks/`)

