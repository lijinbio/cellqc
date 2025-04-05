rule postproc:
  input:
    "result/{sample}.h5ad",
    "nuclear_fraction/{sample}.txt.gz",
  output:
    "postproc/{sample}.h5ad",
  params:
    sampleid="{sample}",
  script:
    "../scripts/postproc.py"
