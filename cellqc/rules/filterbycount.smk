rule filterbycount:
  input:
    "soupx/{sample}.h5" if config["dropkick"]["skip"] else "h5subset/{sample}.h5",
  output:
    directory("filterbycount/{sample}"),
    "filterbycount/{sample}.h5seurat",
    "filterbycount/{sample}/filter_ncell.txt",
  params:
    mincount=config["filterbycount"]["mincount"],
    minfeature=config["filterbycount"]["minfeature"],
    mito=config["filterbycount"]["mito"],
    sampleid="{sample}",
  script:
    "../scripts/filterbycount.R"
