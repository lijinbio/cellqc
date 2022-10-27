rule qcreport:
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
