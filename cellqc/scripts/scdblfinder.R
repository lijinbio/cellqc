# vim: set noexpandtab tabstop=2:
#
# scDblFinder. Diagnostic by default: it scores every cell, but only the caller
# named by doublet.decider removes any. It is fed the SAME expected doublet rate
# as DoubletFinder so the two are compared under matched assumptions rather than
# each using its own default -- otherwise a difference in calls could just be a
# difference in priors.

suppressPackageStartupMessages({
	library(zellkonverter)
	library(SingleCellExperiment)
	library(scDblFinder)
	library(ggplot2)
})

infile=snakemake@input[[1]]
outmeta=snakemake@output[[1]]
outratio=snakemake@output[[2]]
outpdf=snakemake@output[[3]]
outpng=snakemake@output[[4]]

sampleid=snakemake@params[['sampleid']]
nreaction=as.numeric(snakemake@params[['nreaction']])
rate=as.numeric(snakemake@params[['rate']])
capacity=as.numeric(snakemake@params[['capacity']])
seed=as.integer(snakemake@params[['seed']])

set.seed(seed)

sce=readH5AD(infile, reader='R', verbose=FALSE)
# scDblFinder's .checkSCE() hard-requires an assay literally named 'counts';
# zellkonverter calls it 'X'.
if (!'counts' %in% assayNames(sce)) {
	assayNames(sce)[assayNames(sce)=='X']='counts'
}
stopifnot('counts' %in% assayNames(sce))

ncell=ncol(sce)
dbr=rate*ncell/(nreaction*capacity)
cat(sprintf('[scdblfinder] %s: %d genes x %d cells, dbr=%.4f (matched to DoubletFinder)\n',
	sampleid, nrow(sce), ncell, dbr))

res=scDblFinder(sce, dbr=dbr, verbose=FALSE)

out=data.frame(
	barcode=colnames(res),
	scdblfinder_score=as.numeric(res$scDblFinder.score),
	scdblfinder_class=ifelse(as.character(res$scDblFinder.class)=='doublet', 'Doublet', 'Singlet'),
	stringsAsFactors=FALSE
	)
utils::write.table(out, file=gzfile(outmeta), quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)

ndoublet=sum(out$scdblfinder_class=='Doublet')
utils::write.table(
	data.frame(
		sampleid=sampleid, caller='scdblfinder', pK=NA_real_, nreaction=nreaction,
		rate=rate, capacity=capacity, doubletratio=round(dbr, 4),
		ncell_before=ncell, nExp=as.integer(round(dbr*ncell)), ndoublet=ndoublet,
		ncell_after=ncell-ndoublet, homotypic_modelled=FALSE
		),
	file=outratio, quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE
	)

theme_cellqc=theme_bw()+theme(
	plot.background=element_blank(), panel.grid.minor=element_blank(),
	panel.border=element_blank(), plot.title=element_text(hjust=0.5, size=10),
	axis.line=element_line(color='black'), axis.text=element_text(color='black')
	)

pdat=data.frame(score=out$scdblfinder_score,
	class=factor(out$scdblfinder_class, levels=c('Singlet', 'Doublet')))
p=ggplot(pdat, aes(x=class, y=score, fill=class))+
	geom_violin(trim=TRUE, show.legend=FALSE)+
	geom_boxplot(width=0.1, fill='white', outlier.shape=NA)+
	scale_fill_manual(values=c(Singlet='#4c72b0', Doublet='#c44e52'))+
	labs(x=NULL, y='scDblFinder score',
		title=sprintf('%s: scDblFinder score (dbr=%.3f)', sampleid, dbr))+
	theme_cellqc
ggsave(p, file=outpdf, width=4.5, height=4, units='in', device=cairo_pdf)
ggsave(p, file=outpng, width=4.5, height=4, units='in', dpi=300)

cat(sprintf('[scdblfinder] %s: %d/%d cells called doublet (%.2f%%)\n',
	sampleid, ndoublet, ncell, 100*ndoublet/ncell))
