rule soupx:
  input:
    get_cellranger,
  output:
    report(
      "soupx/{sample}_rho.pdf",
      caption="../report/soupx.rst",
      category="Step 1: Ambient RNA removal",
    ),
    report(
      "soupx/{sample}_rhoEst.txt",
      caption="../report/soupx.rst",
      category="Step 1: Ambient RNA removal",
    ),
    "soupx/{sample}.h5",
  params:
    sampleid="{sample}",
  script:
    "../scripts/soupx.R"
