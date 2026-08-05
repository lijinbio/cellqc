rule scdblfinder:
  input:
    "filterbycount/{sample}.h5ad",
  output:
    "scdblfinder/{sample}_metadata.txt.gz",
    "scdblfinder/{sample}_doublet_ratio.txt",
    "scdblfinder/{sample}_score.pdf",
    "scdblfinder/{sample}_score.png",
  params:
    sampleid="{sample}",
    nreaction=get_nreaction,
    rate=config["doublet"]["rate"],
    capacity=config["doublet"]["capacity"],
    seed=config["seed"],
  script:
    "../scripts/scdblfinder.R"
