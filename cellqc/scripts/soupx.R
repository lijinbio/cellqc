suppressPackageStartupMessages(library(SoupX))
suppressPackageStartupMessages(library(DropletUtils))

x=load10X(snakemake@input[[1]])
pdf(snakemake@output[[1]], width=6, height=4)
x=autoEstCont(x)
dev.off()
utils::write.table(data.frame(sampleid=snakemake@params[['sampleid']], rhoEst=x$fit$rhoEst), file=snakemake@output[[2]], quote=F, sep='\t', row.names=F, col.names=T)

x=adjustCounts(x, roundToInt=T)
write10xCounts(snakemake@output[[3]], x)
