rule dropkick:
  input:
    get_rawh5,
  output:
    report(
      "dropkick/{sample}_qc_summary.png",
      caption="../report/dropkick-qc.rst",
      category="Step 2: Empty droplet removal",
    ),
    report(
      "dropkick/{sample}_score_plot.png",
      caption="../report/dropkick-score.rst",
      category="Step 2: Empty droplet removal",
    ),
    report(
      "dropkick/{sample}_coef_plot.png",
      caption="../report/dropkick-coef.rst",
      category="Step 2: Empty droplet removal",
    ),
    "dropkick/{sample}.h5ad",
    "dropkick/{sample}_obs.txt.gz",
    "dropkick/{sample}_var.txt.gz",
  params:
    method=config["dropkick"]["method"],
    numthreads=config["dropkick"]["numthreads"],
  script:
    "../scripts/dropkick.py"
