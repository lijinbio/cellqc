# Terminus of the QC chain. Every caller in doublet.run contributes its
# score/class to .obs; only doublet.decider removes cells. Which callers ran is
# a config value rather than a branch in the DAG, which is why the old three-way
# conditional chaining is gone.
#
# The matrix lands in filterdoublet/ like every other stage's output, and so do
# the doublet statistics; `publish` links those into result/ for the samples the
# qcstatus checkpoint includes. `result/` holds what a user takes away: the final
# matrix from postproc, the doublet statistics, and the reports.
rule filterdoublet:
  input:
    h5ad="filterbycount/{sample}.h5ad",
    metadata=doublet_metadata,
    # Every caller's status, in doublet_callers order: a non-deciding caller
    # that failed is tolerated, so its status is read rather than required.
    caller_status=doublet_status,
    upstream=["filterbycount/{sample}_status.tsv", f"{config['doublet']['decider']}/{{sample}}_status.tsv"],
  output:
    # temp(): postproc is the only consumer, and its output is this matrix plus
    # the sample-prefixed barcodes, unique var names and nuclear fraction. Keeping
    # both meant every run wrote the count matrix to disk twice (72 MB per sample
    # on the reference cohort) for a file whose only unique content was the
    # unprefixed barcode. `snakemake --notemp` keeps it; the cellqc CLI does not
    # pass that through.
    h5ad=temp("filterdoublet/{sample}.h5ad"),
    summary="filterdoublet/{sample}_doublet_summary.txt",
    concordance="filterdoublet/{sample}_doublet_concordance.txt",
    status="filterdoublet/{sample}_status.tsv",
  params:
    sampleid="{sample}",
    callers=doublet_callers,
    decider=config["doublet"]["decider"],
  script:
    "../scripts/filterdoublet.py"
