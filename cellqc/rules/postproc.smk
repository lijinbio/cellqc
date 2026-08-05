rule postproc:
  input:
    unpack(postproc_input),
  output:
    "postproc/{sample}.h5ad",
  params:
    sampleid="{sample}",
  script:
    "../scripts/postproc.py"
