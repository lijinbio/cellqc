rule scpred:
    input:
        "doubletfinder/{sample}.h5",
    output:
        directory("scpred/{sample}"),
        "scpred/{sample}.h5seurat",
    params:
        ref=config["scpred"]["reference"],
        threshold=config["scpred"]["threshold"],
    script:
        "../scripts/scpred.R"
