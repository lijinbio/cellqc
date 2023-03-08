#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source "$(conda info --base)/etc/profile.d/conda.sh"
mamba create -y -n cellqc_v0.0.6 cellqc numpy=1.21 anndata=0.7.8 scanpy=1.9.1 matplotlib=3.6
conda activate cellqc_v0.0.6
Rscript -e "remotes::install_github(c('chris-mcginnis-ucsf/DoubletFinder', 'mojaveazure/seurat-disk'))"
pip install dropkick
pip install git+https://github.com/lijinbio/cellqc.git
