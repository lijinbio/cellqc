rule filterbycount:
  input:
    # Both matrices: the corrected one is what gets filtered and written, the
    # Cell Ranger one only supplies the pre-correction (`raw_*`) QC metrics.
    corrected="ambient/{sample}.h5",
    raw=get_filteredh5,
  output:
    "filterbycount/{sample}.h5ad",
    "filterbycount/{sample}_violin_before.pdf",
    "filterbycount/{sample}_violin_before.png",
    "filterbycount/{sample}_violin_after.pdf",
    "filterbycount/{sample}_violin_after.png",
    "filterbycount/{sample}_filter_ncell.txt",
  params:
    mincount=config["filterbycount"]["mincount"],
    minfeature=config["filterbycount"]["minfeature"],
    mito=config["filterbycount"]["mito"],
    sampleid="{sample}",
    seed=config["seed"],
  script:
    "../scripts/filterbycount.py"
