rule qcreport:
    input:
        soupxrhoEst=expand(["soupx/{sample}_rhoEst.txt"], sample=samples["sample"].tolist()),
        soupxrho=expand(["soupx/{sample}_rho.pdf"], sample=samples["sample"].tolist()),
        filterbycountncell=expand(["filterbycount/{sample}/filter_ncell.txt"], sample=samples["sample"].tolist()),
        filterbycountpltbf=expand(["filterbycount/{sample}/feature_bf.pdf"], sample=samples["sample"].tolist()),
        filterbycountpltaf=expand(["filterbycount/{sample}/feature_af.pdf"], sample=samples["sample"].tolist()),
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
