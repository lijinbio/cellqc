# Final stage: the matrix a user takes away. Written to postproc/ as temp() and
# published to result/ by `publish` (qcstatus.smk) only for the samples the
# qcstatus checkpoint includes, so a failed sample's placeholder never lands in
# result/. The hard link costs no second copy, and the temp() name is removed
# once it is published.
rule postproc:
  input:
    unpack(postproc_input),
  output:
    h5ad=temp("postproc/{sample}.h5ad"),
    obs=temp("postproc/{sample}_obs.txt.gz"),
    var=temp("postproc/{sample}_var.txt.gz"),
    status="postproc/{sample}_status.tsv",
  params:
    sampleid="{sample}",
  script:
    "../scripts/postproc.py"
