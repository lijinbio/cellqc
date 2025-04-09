mamba create -y -n cellqc python=3.10 cellqc r-seurat=4 r-seuratobject=4 r-matrix=1.6.1 dropkick r-hdf5r hdf5 r-leidenbase libxml2 r-xml r-xml2 zlib bioconductor-rsamtools
conda activate cellqc
Rscript -e "remotes::install_github(c('mojaveazure/seurat-disk', 'immunogenomics/harmony', 'powellgenomicslab/scPred', 'powellgenomicslab/DropletQC'), upgrade=F)"
Rscript -e "remotes::install_github('lijinbio/DoubletFinder', upgrade=F, force=T)"
