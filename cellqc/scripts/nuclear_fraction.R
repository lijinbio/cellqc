suppressPackageStartupMessages(library(DropletQC))

crdir=snakemake@input[[1]]
outfile=snakemake@output[[1]]
sampleid=snakemake@params[['sampleid']]
cbtag=snakemake@params[['cbtag']]
retag=snakemake@params[['retag']]
numthreads=snakemake@params[['numthreads']]

result=nuclear_fraction_tags(
	outs=crdir,
	cores=numthreads,
	tiles=100,
	cell_barcode_tag=cbtag,
	region_type_tag=retag,
	verbose=F
)

utils::write.table(
	data.frame(barcode=rownames(result), result),
	file=gzfile(outfile),
	quote=F,
	sep='\t',
	row.names=F,
	col.names=T
	)
