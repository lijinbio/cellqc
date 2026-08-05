rule qcreport:
  input:
    **report_inputs(),
  output:
    html="result/report.html",
    # Machine-readable twin of the report: every scalar, one row per sample.
    metrics="result/metrics.csv",
  params:
    samples=lambda wildcards: samples,
    sampledir=sampledir,
    config=lambda wildcards: config,
    nf_samples=lambda wildcards: nf_samples,
    callers=doublet_callers,
  script:
    "../scripts/qcreport.py"
