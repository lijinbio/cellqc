rule filterbycount:
  input:
    "soupx/{sample}.h5" if config["dropkick"]["skip"] else "h5subset/{sample}.h5",
  output:
    report(
      directory("filterbycount/{sample}"),
      patterns=["{name}.pdf"],
      caption="../report/filterbycount.rst",
      category="Step 3: Filter cells by counts",
    ),
    "filterbycount/{sample}.h5seurat",
    report(
      "filterbycount/{sample}/filter_ncell.txt",
      caption="../report/filterbycount.rst",
      category="Step 3: Filter cells by counts",
    ),
  params:
    mincount=config["filterbycount"]["mincount"],
    minfeature=config["filterbycount"]["minfeature"],
    mito=config["filterbycount"]["mito"],
    sampleid="{sample}",
  script:
    "../scripts/filterbycount.R"
