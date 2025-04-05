rule soupx:
  input:
    get_cellranger,
  output:
    "soupx/{sample}.txt.gz",
  params:
    sampleid="{sample}",
    cbtag=config["nuclear_fraction"]["cbtag"],
    retag=config["nuclear_fraction"]["retag"],
    numthreads=config["nuclear_fraction"]["numthreads"],
  script:
    "../scripts/nuclear_fraction.R"
