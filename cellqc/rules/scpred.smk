rule scpred:
    input:
        "doubletfinder/{sample}.rds",
    output:
        report(
            directory("scpred/{sample}"),
            patterns=["{name1}.png", "{name2}.txt"],
            caption="../report/scpred.rst",
            category="Step 4: Cell type annotation",
        ),
        "scpred/{sample}.h5seurat",
    params:
        ref=config["scpred"]["reference"],
        threshold=config["scpred"]["threshold"],
    script:
        "../scripts/scpred.R"
