rule soupx:
  input:
    get_cellranger,
  output:
    "soupx/{sample}_rho.pdf",
    "soupx/{sample}_rhoEst.txt",
    "soupx/{sample}.h5",
  params:
    sampleid="{sample}",
  script:
    "../scripts/soupx.R"
