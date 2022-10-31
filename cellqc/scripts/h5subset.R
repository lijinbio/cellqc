suppressPackageStartupMessages(library(DropletUtils))

x=read10xCounts(snakemake@input[[1]], col.names=T)
f=read.table(snakemake@input[[2]], header=T)
f=subset(f, dropkick_label=='True')
ids=f[, 'barcode']

res=data.frame(
	sampleid=snakemake@params[['sampleid']]
	, emptydrops=ncol(x)
	, dropkick=length(ids)
	, intersection=length(intersect(ids, colnames(x)))
	)
write.table(res, file=snakemake@output[[1]], quote=F, sep='\t', row.names=F, col.names=T)

x=x[, intersect(ids, colnames(x)), drop=F]
x=counts(x)
write10xCounts(snakemake@output[[2]], x)
