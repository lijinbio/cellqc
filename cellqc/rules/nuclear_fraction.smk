# Enabled per sample by the presence of an indexed BAM (see rules/config.smk).
# There is no skip flag: `final_targets()` only ever asks for the samples in
# `nf_samples`, so a cohort with mixed BAM availability just works.
rule nuclear_fraction:
  input:
    cellranger=get_cellranger,
    filtered=get_filteredh5,
  output:
    "nuclear_fraction/{sample}.txt.gz",
    "nuclear_fraction/{sample}_nf_umi.pdf",
    "nuclear_fraction/{sample}_nf_umi.png",
  params:
    sampleid="{sample}",
    cbtag=config["nuclear_fraction"]["cbtag"],
    retag=config["nuclear_fraction"]["retag"],
    exontag=config["nuclear_fraction"]["exontag"],
    introntag=config["nuclear_fraction"]["introntag"],
  threads:
    config["nuclear_fraction"]["numthreads"]
  script:
    "../scripts/nuclear_fraction.py"
