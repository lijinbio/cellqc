# Final stage: the matrix a user takes away. Written to postproc/ as temp() and
# published to result/ by `publish` (qcstatus.smk) only for the samples the
# qcstatus checkpoint includes, so a failed sample's placeholder never lands in
# result/. The hard link costs no second copy, and the temp() name is removed
# once it is published. The status file is temp() too: the qcstatus checkpoint
# is its only reader and copies it into result/qc_status.csv, so nothing is left
# in postproc/ and Snakemake removes the empty directory at the end of the run.
rule postproc:
  input:
    unpack(postproc_input),
  output:
    h5ad=temp("postproc/{sample}.h5ad"),
    obs=temp("postproc/{sample}_obs.txt.gz"),
    var=temp("postproc/{sample}_var.txt.gz"),
    status=temp("postproc/{sample}_status.tsv"),
  params:
    sampleid="{sample}",
  script:
    "../scripts/postproc.py"
