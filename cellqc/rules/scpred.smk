rule scpred:
  input:
    "doubletfinder/{sample}.h5seurat" if not config['doubletfinder']['skip'] else "filterbycount/{sample}.h5seurat",
  output:
    directory("scpred/{sample}"),
    "result/{sample}.h5seurat",
    "scpred/{sample}/contingency.txt",
  params:
    ref=config["scpred"]["reference"],
    threshold=config["scpred"]["threshold"],
  script:
    "../scripts/scpred.R"
