rule h5subset:
  input:
    "soupx/{sample}.h5",
    "dropkick/{sample}_obs.txt.gz",
  output:
    report(
      "h5subset/{sample}_dropkick_stat.txt",
      caption="../report/h5subset.rst",
      category="Step 2: Empty droplet removal",
    ),
    "h5subset/{sample}.h5",
  params:
    sampleid="{sample}",
  script:
    "../scripts/h5subset.R"
