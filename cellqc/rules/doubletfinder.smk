rule doubletfinder:
    input:
        "filterbycount/{sample}.h5seurat",
    output:
        report(
            directory("doubletfinder/{sample}"),
            patterns=["{name}.pdf"],
            caption="../report/doubletfinder.rst",
            category="Step 4: Doublet removal",
        ),
        "doubletfinder/{sample}.rds",
    params:
        findpK=config["doubletfinder"]["findpK"],
        numthreads=config["doubletfinder"]["numthreads"],
        pK=config["doubletfinder"]["pK"],
        nrun=lambda wildcards: samples.loc[wildcards.sample, "nrun"]
        sampleid="{sample}"
    script:
        "../scripts/doubletfinder.R"
