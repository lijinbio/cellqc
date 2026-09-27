rule ambient:
  input:
    get_cellranger,
  output:
    "ambient/{sample}.h5",
    "ambient/{sample}_contamination.txt",
    "ambient/{sample}_ambient.pdf",
    "ambient/{sample}_ambient.png",
    status="ambient/{sample}_status.tsv",
  params:
    sampleid="{sample}",
    method=config["ambient"]["method"],
    compare=config["ambient"]["compare"],
    seed=config["seed"],
    guard=GUARD_R,
  script:
    "../scripts/ambient.R"
