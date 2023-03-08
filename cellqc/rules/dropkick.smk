rule dropkick:
  input:
    get_rawh5,
  output:
    "dropkick/{sample}_qc_summary.png",
    "dropkick/{sample}.h5ad",
  params:
    method=config["dropkick"]["method"],
    numthreads=config["dropkick"]["numthreads"],
  script:
    "../scripts/dropkick.py"
