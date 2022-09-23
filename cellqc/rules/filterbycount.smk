rule filterbycount:
    input:
        "scpred/{sample}.h5seurat",
    output:
        report(
            directory("filterbycount/{sample}"),
            patterns=["{name1}.png", "{name2}.txt"],
            caption="../report/filterbycount.rst",
            category="Filter cells by counts",
        ),
        "filterbycount/{sample}.h5seurat",
    params:
        mincount=config["filterbycount"]["mincount"],
        minfeature=config["filterbycount"]["minfeature"],
        mito=config["filterbycount"]["mito"],
        sample="{sample}"
    script:
        "../scripts/filterbycount.R"
