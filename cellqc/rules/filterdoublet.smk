# Single terminus of the QC chain. Every caller in doublet.run contributes its
# score/class to .obs; only doublet.decider removes cells. When doublet.skip is
# set there is no metadata at all and this rule just relabels the matrix, which
# is why the old three-way conditional chaining is gone.
rule filterdoublet:
  input:
    h5ad="filterbycount/{sample}.h5ad",
    metadata=doublet_metadata,
  output:
    "result/{sample}.h5ad",
    "result/{sample}_doublet_summary.txt",
    "result/{sample}_doublet_concordance.txt",
  params:
    sampleid="{sample}",
    callers=doublet_callers,
    decider=config["doublet"]["decider"],
    skip=config["doublet"]["skip"],
  script:
    "../scripts/filterdoublet.py"
