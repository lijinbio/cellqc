rule soupx:
    input:
        get_cellranger,
    output:
        report(
            "soupx/{sample}_rho.pdf",
            caption="../report/soupx.rst",
            category="Ambient RNA removal",
        ),
        "soupx/{sample}_rhoEst.txt",
        "soupx/{sample}.h5",
    script:
        "../scripts/soupx.R"
