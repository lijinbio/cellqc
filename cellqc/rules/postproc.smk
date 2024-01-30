rule postproc:
  input:
    "result/{sample}.h5ad",
  output:
    "postproc/{sample}.h5ad",
  params:
    sampleid="{sample}",
  script:
    "../scripts/postproc.py"
