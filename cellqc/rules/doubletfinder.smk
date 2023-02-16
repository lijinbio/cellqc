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
    "result/{sample}.h5seurat" if config["scpred"]["skip"] else "doubletfinder/{sample}.h5seurat",
    report(
      "doubletfinder/{sample}/doublet_ratio.txt",
      caption="../report/doubletfinder.rst",
      category="Step 4: Doublet removal",
    ),
  params:
    findpK=config["doubletfinder"]["findpK"],
    numthreads=config["doubletfinder"]["numthreads"],
    pK=config["doubletfinder"]["pK"],
    nreaction=lambda wildcards: samples.loc[wildcards.sample, "nreaction"] if "nreaction" in samples.columns else 1,
    sampleid="{sample}",
  script:
    "../scripts/doubletfinder.R"
