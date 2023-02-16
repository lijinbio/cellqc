rule qcreport:
  input:
    soupxrhoEst=expand(["soupx/{sample}_rhoEst.txt"], sample=samples["sample"].tolist()),
    soupxrho=expand(["soupx/{sample}_rho.pdf"], sample=samples["sample"].tolist()),
    dropkickstat=expand(["h5subset/{sample}_dropkick_stat.txt"], sample=samples["sample"].tolist() if not config["dropkick"]["skip"] else []),
    dropkicksummary=expand(["dropkick/{sample}_qc_summary.png"], sample=samples["sample"].tolist() if not config["dropkick"]["skip"] else []),
    dropkickscore=expand(["dropkick/{sample}_score_plot.png"], sample=samples["sample"].tolist() if not config["dropkick"]["skip"] else []),
    filterbycountdir=expand(["filterbycount/{sample}"], sample=samples["sample"].tolist()),
    doubletfinderdir=expand(["doubletfinder/{sample}"], sample=samples["sample"].tolist()),
    scpreddir=expand(["scpred/{sample}"], sample=samples["sample"].tolist() if not config["scpred"]["skip"] else []),
  output:
    report(
      "result/qc_report.html",
      caption="../report/qcreport.rst",
      category="Step 6: QC report",
    ),
  params:
    samples=samples,
    sampledir=sampledir,
    pK=config["doubletfinder"]["pK"],
  script:
    "../scripts/qcreport.py"
