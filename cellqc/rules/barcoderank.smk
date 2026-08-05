rule barcoderank:
  input:
    raw=get_rawh5,
    filtered=get_filteredh5,
  output:
    "barcoderank/{sample}_barcoderank.pdf",
    "barcoderank/{sample}_barcoderank.png",
    "barcoderank/{sample}_knee.txt",
  params:
    sampleid="{sample}",
  script:
    "../scripts/barcoderank.py"
