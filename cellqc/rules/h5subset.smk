rule h5subset:
  input:
    "soupx/{sample}.h5",
    "dropkick/{sample}_obs.txt.gz",
  output:
    "h5subset/{sample}_dropkick_stat.txt",
    "h5subset/{sample}.h5",
  params:
    sampleid="{sample}",
  script:
    "../scripts/h5subset.R"
