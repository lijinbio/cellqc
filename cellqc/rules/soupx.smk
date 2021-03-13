rule soupx:
    input:
        get_cellranger,
    output:
        "soupx/{sample}_rho.pdf",
        "soupx/{sample}.h5",
    script:
        "../scripts/soupx.R"
