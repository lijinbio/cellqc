# Final stage: the matrix a user takes away, in result/.
rule postproc:
  input:
    unpack(postproc_input),
  output:
    h5ad="result/{sample}.h5ad",
    obs="result/{sample}_obs.txt.gz",
    var="result/{sample}_var.txt.gz",
  params:
    sampleid="{sample}",
  script:
    "../scripts/postproc.py"
