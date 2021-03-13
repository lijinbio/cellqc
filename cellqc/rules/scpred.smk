rule scpred:
    input:
        "doubletfinder/{sample}.h5",
    output:
        report(
            directory("scpred/{sample}"),
            patterns=["*.png", "*.txt"],
            caption="../report/scpred.rst",
            category="Cell type annotation",
        ),
        "scpred/{sample}.h5seurat",
    params:
        ref=config["scpred"]["reference"],
        threshold=config["scpred"]["threshold"],
    script:
        "../scripts/scpred.R"
