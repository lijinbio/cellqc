rule h5seurat2h5ad:
  input:
    "result/{sample}.h5seurat",
  output:
    "result/{sample}.h5ad",
  script:
    "../scripts/h5seurat2h5ad.R"
