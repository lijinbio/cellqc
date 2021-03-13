rule h5seurat2h5ad:
    input:
        "scpred/{sample}.h5seurat",
    output:
        "cellqc/{sample}.h5ad",
    script:
        "../scripts/h5seurat2h5ad.R"
