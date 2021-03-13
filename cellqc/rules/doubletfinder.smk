rule doubletfinder:
    input:
        "h5subset/{sample}.h5",
    output:
        report(
            directory("doubletfinder/{sample}"),
            patterns=["*.pdf"],
            caption="../report/doubletfinder.rst",
            category="Doublet removal",
        ),
        "doubletfinder/{sample}.h5",
    params:
        findpK=config["doubletfinder"]["findpK"],
        numthreads=config["doubletfinder"]["numthreads"],
        pK=config["doubletfinder"]["pK"],
    script:
        "../scripts/doubletfinder.R"
