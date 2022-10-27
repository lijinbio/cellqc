rule scpred:
    input:
        "doubletfinder/{sample}.h5seurat" if not config["scpred"]["skip"],
    output:
        report(
            directory("scpred/{sample}"),
            patterns=["{name}.png"],
            caption="../report/scpred.rst",
            category="Step 5: Cell type annotation",
        ),
        "result/{sample}.h5seurat",
        report(
            "scpred/{sample}/contingency.txt",
            caption="../report/scpred.rst",
            category="Step 5: Cell type annotation",
        ),
    params:
        ref=config["scpred"]["reference"],
        threshold=config["scpred"]["threshold"],
    script:
        "../scripts/scpred.R"
