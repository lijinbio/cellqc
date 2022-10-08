rule qcreport:
    output:
        report(
            "qc_report.html",
            caption="../report/qcreport.rst",
            category="Step 6: QC report",
        ),
    script:
        "../scripts/qcreport.py"
