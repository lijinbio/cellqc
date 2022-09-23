rule doubletfinder:
    input:
        "h5subset/{sample}.h5",
    output:
        report(
            directory("doubletfinder/{sample}"),
            patterns=["{name}.pdf"],
            caption="../report/doubletfinder.rst",
            category="Step 3: Doublet removal",
        ),
        "doubletfinder/{sample}.rds",
    params:
        findpK=config["doubletfinder"]["findpK"],
        numthreads=config["doubletfinder"]["numthreads"],
        pK=config["doubletfinder"]["pK"],
        nrun=lambda wildcards: samples.loc[wildcards.sample, "nrun"]
    script:
        "../scripts/doubletfinder.R"
