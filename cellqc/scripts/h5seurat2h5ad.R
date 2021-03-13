# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(SeuratDisk))
Convert(snakemake@input[[1]], dest=snakemake@output[[1]])
