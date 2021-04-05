#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source "$(conda info --base)/etc/profile.d/conda.sh"
conda init && conda create -y -n cellqc && conda activate cellqc
conda config --add channels defaults --add channels bioconda --add channels conda-forge
set -e
conda install -y mamba
mamba install -y bioconductor-dropletutils r-seurat r-dplyr r-ggplot2 r-soupx r-remotes scanpy pygraphviz snakemake
Rscript -e "remotes::install_github(c('chris-mcginnis-ucsf/DoubletFinder', 'mojaveazure/seurat-disk', 'immunogenomics/harmony', 'powellgenomicslab/scPred'))"
pip install dropkick
pip install .
