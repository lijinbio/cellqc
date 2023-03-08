rule doubletfinder:
  input:
    "filterbycount/{sample}.h5seurat",
  output:
    directory("doubletfinder/{sample}"),
    "doubletfinder/{sample}/doublet_ratio.txt",
    "result/{sample}.h5seurat" if config["scpred"]["skip"] else "doubletfinder/{sample}.h5seurat",
  params:
    findpK=config["doubletfinder"]["findpK"],
    numthreads=config["doubletfinder"]["numthreads"],
    pK=config["doubletfinder"]["pK"],
    nreaction=lambda wildcards: samples.loc[wildcards.sample, "nreaction"] if "nreaction" in samples.columns else 1,
    sampleid="{sample}",
  script:
    "../scripts/doubletfinder.R"
