#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source "$(conda info --base)/etc/profile.d/conda.sh"
conda create -y -n cellqc && conda activate cellqc
conda config --add channels defaults --add channels bioconda --add channels conda-forge
set -ex
conda install -y mamba
mamba install -y bioconductor-dropletutils r-seurat r-seuratobject r-dplyr r-ggplot2 'r-soupx>=1.6.2' r-remotes scanpy pygraphviz snakemake click
Rscript -e "remotes::install_github(c('chris-mcginnis-ucsf/DoubletFinder', 'mojaveazure/seurat-disk', 'immunogenomics/harmony', 'powellgenomicslab/scPred'))"
mamba install -y numpy=1.21 # bug fix to install dropkick
mamba install -y anndata=0.7.8 # Fix .h5ad version
pip install dropkick
pip install git+https://github.com/lijinbio/cellqc.git
