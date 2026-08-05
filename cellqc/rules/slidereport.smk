rule slidereport:
  input:
    **report_inputs(),
  output:
    "result/report_slides.pdf",
  params:
    samples=lambda wildcards: samples,
    sampledir=sampledir,
    config=lambda wildcards: config,
    nf_samples=lambda wildcards: nf_samples,
    callers=doublet_callers,
    nowtimestr=nowtimestr,
  script:
    "../scripts/slidereport.py"
