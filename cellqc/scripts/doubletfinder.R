# vim: set noexpandtab tabstop=2:
#
# DoubletFinder. Reads the .h5ad written by filterbycount.py through
# zellkonverter's native R reader (no reticulate/basilisk), and writes only a
# per-barcode metadata table -- filterdoublet.py applies it. R never rewrites
# the matrix, which keeps a second serializer out of the count path.

suppressPackageStartupMessages({
	library(zellkonverter)
	library(SingleCellExperiment)
	library(Seurat)
	library(DoubletFinder)
	library(ggplot2)
})

infile=snakemake@input[[1]]
outmeta=snakemake@output[[1]]
outratio=snakemake@output[[2]]
outpann_pdf=snakemake@output[[3]]
outpann_png=snakemake@output[[4]]
outumap_pdf=snakemake@output[[5]]
outumap_png=snakemake@output[[6]]

sampleid=snakemake@params[['sampleid']]
findpK=snakemake@params[['findpK']]
pK=snakemake@params[['pK']]
nreaction=as.numeric(snakemake@params[['nreaction']])
rate=as.numeric(snakemake@params[['rate']])
capacity=as.numeric(snakemake@params[['capacity']])
seed=as.integer(snakemake@params[['seed']])
numthreads=as.integer(snakemake@threads)

# v0.1.0 seeded nothing: RunPCA/RunUMAP and DoubletFinder's artificial-doublet
# sampling are all stochastic, so repeat runs gave different calls.
set.seed(seed)

sce=readH5AD(infile, reader='R', verbose=FALSE)
# zellkonverter names the main assay 'X'; Seurat and scDblFinder both expect
# 'counts'. Normalise the name once, here.
if (!'counts' %in% assayNames(sce)) {
	assayNames(sce)[assayNames(sce)=='X']='counts'
}
stopifnot('counts' %in% assayNames(sce))

x=CreateSeuratObject(counts=assay(sce, 'counts'), meta.data=as.data.frame(colData(sce)))
cat(sprintf('[doubletfinder] %s: %d genes x %d cells\n', sampleid, nrow(x), ncol(x)))

npc=min(30, ncol(x)-1)
x=NormalizeData(x, verbose=FALSE)
x=FindVariableFeatures(x, selection.method='vst', nfeatures=2000, verbose=FALSE)
x=ScaleData(x, verbose=FALSE)
x=RunPCA(x, npcs=npc, verbose=FALSE, seed.use=seed)
x=RunUMAP(x, dims=1:npc, verbose=FALSE, seed.use=seed)

## ---- pK --------------------------------------------------------------------
if (isTRUE(findpK)) {
	sweepx=paramSweep(x, PCs=1:10, sct=FALSE, num.cores=numthreads)
	sweepstats=summarizeSweep(sweepx, GT=FALSE)
	bcmetric=find.pK(sweepstats)
	utils::write.table(bcmetric, file=gzfile(sub('_metadata.txt.gz$', '_findpK.txt.gz', outmeta)),
		quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
	bcmetric$pK=as.numeric(levels(bcmetric$pK))[bcmetric$pK]
	pKopt=bcmetric$pK[which.max(bcmetric[, 'BCmetric'])]
	cat(sprintf('[doubletfinder] %s: estimated pK=%g\n', sampleid, pKopt))
} else {
	pKopt=pK
	cat(sprintf('[doubletfinder] %s: using preset pK=%g\n', sampleid, pKopt))
}

## ---- Expected doublets -----------------------------------------------------
# 10x multiplet rule of thumb, unchanged from v0.1.0 but with the two constants
# now configurable. NOT adjusted for homotypic doublets (modelHomotypic is
# deliberately not called), so nExp over-estimates the DETECTABLE doublet count
# and this step removes somewhat more cells than the true heterotypic count.
# The bias direction is known, constant, and reported.
ncell=ncol(x)
doubletratio=round(rate*ncell/(nreaction*capacity), digits=2)
nExp=as.integer(round(doubletratio*ncell))
cat(sprintf('[doubletfinder] %s: expected doublet rate %.3f of %d cells -> nExp=%d (homotypic NOT modelled)\n',
	sampleid, doubletratio, ncell, nExp))

# reuse.pANN must be NULL, not FALSE. Upstream v2.0.6 switched this check from
# `if (reuse.pANN)` to `if (!is.null(reuse.pANN))`, so the FALSE that v0.1.0
# passed now takes the reuse branch and fails with "cannot xtfrm data frames".
res=doubletFinder(x, PCs=1:10, pN=0.25, pK=pKopt, nExp=nExp, reuse.pANN=NULL, sct=FALSE)

meta=res@meta.data
scorecol=grep('^pANN_', names(meta), value=TRUE)[1]
classcol=grep('^DF.classifications_', names(meta), value=TRUE)[1]
stopifnot(!is.na(scorecol), !is.na(classcol))

out=data.frame(
	barcode=rownames(meta),
	doubletfinder_pANN=meta[[scorecol]],
	doubletfinder_class=as.character(meta[[classcol]]),
	stringsAsFactors=FALSE
	)
utils::write.table(out, file=gzfile(outmeta), quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)

ndoublet=sum(out$doubletfinder_class=='Doublet')
utils::write.table(
	data.frame(
		sampleid=sampleid, caller='doubletfinder', pK=pKopt, nreaction=nreaction,
		rate=rate, capacity=capacity, doubletratio=doubletratio,
		ncell_before=ncell, nExp=nExp, ndoublet=ndoublet, ncell_after=ncell-ndoublet,
		homotypic_modelled=FALSE
		),
	file=outratio, quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE
	)

## ---- Figures ---------------------------------------------------------------
res$pANN=out$doubletfinder_pANN
res$DF_class=factor(out$doubletfinder_class, levels=c('Singlet', 'Doublet'))

theme_cellqc=theme_bw()+theme(
	plot.background=element_blank(), panel.grid.minor=element_blank(),
	panel.border=element_blank(), plot.title=element_text(hjust=0.5, size=10),
	axis.line=element_line(color='black'), axis.text=element_text(color='black')
	)

pdat=data.frame(pANN=res$pANN, class=res$DF_class)
p=ggplot(pdat, aes(x=class, y=pANN, fill=class))+
	geom_violin(trim=TRUE, show.legend=FALSE)+
	geom_boxplot(width=0.1, fill='white', outlier.shape=NA)+
	scale_fill_manual(values=c(Singlet='#4c72b0', Doublet='#c44e52'))+
	labs(x=NULL, y='pANN', title=sprintf('%s: pANN by classification (pK=%g)', sampleid, pKopt))+
	theme_cellqc
ggsave(p, file=outpann_pdf, width=4.5, height=4, units='in', device=cairo_pdf)
ggsave(p, file=outpann_png, width=4.5, height=4, units='in', dpi=300)

p=DimPlot(res, reduction='umap', group.by='DF_class', cols=c(Singlet='#4c72b0', Doublet='#c44e52'))+
	labs(title=sprintf('%s: DoubletFinder calls', sampleid))+theme_cellqc
ggsave(p, file=outumap_pdf, width=5, height=4.2, units='in', device=cairo_pdf)
ggsave(p, file=outumap_png, width=5, height=4.2, units='in', dpi=300)

cat(sprintf('[doubletfinder] %s: %d/%d cells called doublet (%.2f%%)\n',
	sampleid, ndoublet, ncell, 100*ndoublet/ncell))
