# Terminus of the QC chain. Every caller in doublet.run contributes its
# score/class to .obs; only doublet.decider removes cells. Which callers ran is
# a config value rather than a branch in the DAG, which is why the old three-way
# conditional chaining is gone.
#
# The matrix lands in filterdoublet/ like every other stage's output. `result/`
# holds what a user takes away: the final matrix from postproc, the doublet
# statistics, and the two reports.
rule filterdoublet:
  input:
    h5ad="filterbycount/{sample}.h5ad",
    metadata=doublet_metadata,
  output:
    h5ad="filterdoublet/{sample}.h5ad",
    obs="filterdoublet/{sample}_obs.txt.gz",
    var="filterdoublet/{sample}_var.txt.gz",
    summary="result/{sample}_doublet_summary.txt",
    concordance="result/{sample}_doublet_concordance.txt",
  params:
    sampleid="{sample}",
    callers=doublet_callers,
    decider=config["doublet"]["decider"],
  script:
    "../scripts/filterdoublet.py"
