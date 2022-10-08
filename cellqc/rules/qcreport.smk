rule qcreport:
    input:
        cellrangersummary=get_cellrangersummary,
    output:
        report(
            "qc_report.html",
            caption="../report/qcreport.rst",
            category="Step 6: QC report",
        ),
    script:
        "../scripts/qcreport.py"
