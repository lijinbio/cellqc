rule qcreport:
    input:
        soupxrhoEst=expand(["soupx/{sample}_rhoEst.txt"], sample=samples["sample"].tolist()),
        soupxrho=expand(["soupx/{sample}_rho.pdf"], sample=samples["sample"].tolist()),
        filterbycountncell=expand(["filterbycount/{sample}/filter_ncell.txt"], sample=samples["sample"].tolist()),
        filterbycountpltbf=expand(["filterbycount/{sample}/feature_bf.pdf"], sample=samples["sample"].tolist()),
        filterbycountpltaf=expand(["filterbycount/{sample}/feature_af.pdf"], sample=samples["sample"].tolist()),
        doubletratio=expand(["doubletfinder/{sample}/doublet_ratio.txt"], sample=samples["sample"].tolist()),
        doubletpannviolin=expand(["doubletfinder/{sample}/pANN_violin_pK{pK}.pdf"], sample=samples["sample"].tolist(), pK=config["doubletfinder"]["pK"]),
        doublettsne=expand(["doubletfinder/{sample}/tsne_doublet_pK{pK}.pdf"], sample=samples["sample"].tolist(), pK=config["doubletfinder"]["pK"]),
        doubletumap=expand(["doubletfinder/{sample}/umap_doublet_pK{pK}.pdf"], sample=samples["sample"].tolist(), pK=config["doubletfinder"]["pK"]),
    output:
        report(
            "result/qc_report.html",
            caption="../report/qcreport.rst",
            category="Step 6: QC report",
        ),
    params:
        samples=samples,
        sampledir=sampledir,
    script:
        "../scripts/qcreport.py"
