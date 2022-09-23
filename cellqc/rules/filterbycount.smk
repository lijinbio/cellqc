rule filterbycount:
    input:
        "scpred/{sample}.h5seurat",
    output:
        report(
            directory("filterbycount/{sample}"),
            patterns=["{name}.pdf"],
            caption="../report/filterbycount.rst",
            category="Step 5: Filter cells by counts",
        ),
        "filterbycount/{sample}.h5seurat",
        report(
            "filterbycount/{sample}/ncell.txt",
            caption="../report/filterbycount.rst",
            category="Step 5: Filter cells by counts",
        ),
    params:
        mincount=config["filterbycount"]["mincount"],
        minfeature=config["filterbycount"]["minfeature"],
        mito=config["filterbycount"]["mito"],
        sample="{sample}"
    script:
        "../scripts/filterbycount.R"
