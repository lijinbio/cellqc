# vim: set noexpandtab tabstop=2:
#
# Ambient RNA estimation and correction.
#
# `method` is the ONE method whose corrected counts are written. Methods in
# `compare` are run for their contamination estimate only and never touch the
# output matrix -- they exist so that disagreement between methods is visible,
# not so a correction can be picked after seeing which one flatters the
# downstream result.

suppressPackageStartupMessages({
	library(SoupX)
	library(DropletUtils)
	library(Matrix)
})

crdir=snakemake@input[[1]]
outh5=snakemake@output[[1]]
outcontam=snakemake@output[[2]]
outpdf=snakemake@output[[3]]
outpng=snakemake@output[[4]]
sampleid=snakemake@params[['sampleid']]
method=snakemake@params[['method']]
compare=snakemake@params[['compare']]
seed=snakemake@params[['seed']]

set.seed(seed)

filt_h5=file.path(crdir, 'filtered_feature_bc_matrix.h5')
raw_h5=file.path(crdir, 'raw_feature_bc_matrix.h5')

## ---- Gene identifiers ------------------------------------------------------
# SoupX::load10X goes through Seurat::Read10X, whose rownames are make.unique()d
# gene SYMBOLS -- the Ensembl IDs are lost. v0.1.0 wrote the matrix that way, so
# every downstream .h5ad carried symbols only. Recover the IDs by rebuilding the
# exact same make.unique() mapping from the feature table, and assert it matches
# rather than trusting it.
feature_table=function() {
	sce=read10xCounts(filt_h5, col.names=TRUE)
	rd=as.data.frame(SummarizedExperiment::rowData(sce))
	idcol=if ('ID' %in% names(rd)) 'ID' else names(rd)[1]
	symcol=if ('Symbol' %in% names(rd)) 'Symbol' else names(rd)[2]
	data.frame(
		gene_id=as.character(rd[[idcol]]),
		gene_symbol=as.character(rd[[symcol]]),
		seurat_rowname=make.unique(as.character(rd[[symcol]])),
		stringsAsFactors=FALSE
		)
}
features=feature_table()

align_features=function(counts, keyed_by) {
	# Return counts reordered to the canonical feature order, plus id/symbol
	# vectors, so every method writes an identically-indexed matrix.
	key=switch(keyed_by, symbol=features$seurat_rowname, id=features$gene_id)
	if (!setequal(rownames(counts), key)) {
		stop(sprintf(
			'Feature mismatch for %s: matrix has %d rows, feature table %d, %d shared. Refusing to write a mislabelled matrix.',
			sampleid, nrow(counts), length(key), length(intersect(rownames(counts), key))))
	}
	counts[key, , drop=FALSE]
}

## ---- Figure ----------------------------------------------------------------
# Drawn from sc$fit rather than letting autoEstCont plot, so the estimate runs
# once and the PDF and PNG are the same figure (autoEstCont is called with
# doPlot=FALSE).
plot_soupx=function(fit) {
	rhoProbes=seq(0, 1, 0.001)
	v2=(fit$priorRhoStdDev/fit$priorRho)^2
	k=1+v2^-2/2*(1+sqrt(1+4*v2))
	theta=fit$priorRho/(k-1)
	prior=dgamma(rhoProbes, k, scale=theta)
	plot(rhoProbes, fit$posterior, type='l', xlim=c(0, 1),
		ylim=c(0, max(c(prior, fit$posterior))), frame.plot=FALSE,
		xlab='Contamination fraction', ylab='Probability density',
		main=sprintf('%s: SoupX rho = %.3f', sampleid, fit$rhoEst))
	lines(rhoProbes, prior, lty=2)
	abline(v=fit$rhoEst, col='red')
	legend('topright', bty='n',
		legend=c(sprintf('prior %g (+/-%g)', fit$priorRho, fit$priorRhoStdDev),
			sprintf('posterior %.3f (%.3f, %.3f)', fit$rhoEst, fit$rhoFWHM[1], fit$rhoFWHM[2])),
		lty=c(2, 1), col=c('black', 'black'))
}

plot_decontx=function(contamination) {
	hist(contamination, breaks=50, col='grey70', border='white',
		main=sprintf('%s: DecontX contamination (mean %.3f)', sampleid, mean(contamination)),
		xlab='Per-cell contamination fraction')
	abline(v=mean(contamination), col='red')
}

draw=function(fn) {
	pdf(outpdf, width=6, height=4.5); fn(); dev.off()
	png(outpng, width=6, height=4.5, units='in', res=200); fn(); dev.off()
}

## ---- SoupX -----------------------------------------------------------------
run_soupx=function() {
	sc=load10X(crdir, verbose=FALSE)
	sc=autoEstCont(sc, doPlot=FALSE, forceAccept=TRUE, verbose=TRUE)
	rho=sc$fit$rhoEst
	list(
		contamination=setNames(rep(rho, ncol(sc$toc)), colnames(sc$toc)),
		global=rho,
		counts=align_features(adjustCounts(sc, roundToInt=TRUE), 'symbol'),
		plot=function() plot_soupx(sc$fit)
		)
}

## ---- DecontX ---------------------------------------------------------------
run_decontx=function() {
	suppressPackageStartupMessages({
		library(celda)
		library(SingleCellExperiment)
	})
	x=read10xCounts(filt_h5, col.names=TRUE)
	stopifnot('counts' %in% assayNames(x))
	bg=if (file.exists(raw_h5)) read10xCounts(raw_h5, col.names=TRUE) else NULL
	res=if (is.null(bg)) decontX(x, seed=seed) else decontX(x, background=bg, seed=seed)
	cont=res$decontX_contamination
	names(cont)=colnames(res)
	list(
		contamination=cont,
		global=mean(cont),
		# decontX returns fractional counts; round so the matrix stays integer
		# UMI counts, which every downstream count-based model assumes.
		counts=align_features(round(assay(res, 'decontXcounts')), 'id'),
		plot=function() plot_decontx(cont)
		)
}

runners=list(soupx=run_soupx, decontx=run_decontx)

methods=unique(c(method, compare))
methods=methods[methods!='none']

estimates=list()
corrected=NULL
for (m in methods) {
	cat(sprintf('[ambient] %s: running %s (%s)\n', sampleid, m,
		if (m==method) 'APPLIED to counts' else 'comparison only, not applied'))
	# Re-seed before EACH method. SoupX::adjustCounts(roundToInt=TRUE) does
	# randomised rounding (rbinom on the fractional part), so the corrected count
	# matrix is stochastic -- v0.1.0 seeded nothing and therefore produced a
	# slightly different integer matrix on every run, which propagated to
	# pct_counts_mt and flipped cells sitting near the mitochondrial threshold.
	# Seeding per method also makes each method's result independent of which
	# other methods run and in what order.
	set.seed(seed)
	estimates[[m]]=runners[[m]]()
	if (m==method) corrected=estimates[[m]]$counts
}

if (method=='none') {
	corrected=align_features(counts(read10xCounts(filt_h5, col.names=TRUE)), 'id')
	draw(function() { plot.new(); title(main=sprintf('%s: no ambient correction applied', sampleid)) })
} else {
	draw(estimates[[method]]$plot)
}

## ---- Write the corrected matrix --------------------------------------------
orig=read10xCounts(filt_h5, col.names=TRUE)
tot_before=sum(counts(orig))
tot_after=sum(corrected)

write10xCounts(
	outh5, corrected,
	barcodes=colnames(corrected),
	gene.id=features$gene_id,
	gene.symbol=features$gene_symbol,
	type='HDF5', version='3', overwrite=TRUE
	)

## ---- Contamination table ---------------------------------------------------
# Long format, one row per method, so applied and comparison-only estimates sit
# side by side and can never be confused for one another.
rows=lapply(names(estimates), function(m) {
	e=estimates[[m]]
	data.frame(
		sampleid=sampleid, method=m, applied=(m==method),
		contamination_mean=mean(e$contamination),
		contamination_median=stats::median(e$contamination),
		contamination_min=min(e$contamination),
		contamination_max=max(e$contamination),
		ncell=length(e$contamination),
		counts_before=tot_before,
		counts_after=ifelse(m==method, tot_after, NA_real_),
		counts_removed_frac=ifelse(m==method, (tot_before-tot_after)/tot_before, NA_real_),
		stringsAsFactors=FALSE
		)
	})
tab=if (length(rows)) do.call(rbind, rows) else data.frame(
	sampleid=sampleid, method='none', applied=TRUE,
	contamination_mean=NA_real_, contamination_median=NA_real_,
	contamination_min=NA_real_, contamination_max=NA_real_,
	ncell=ncol(corrected), counts_before=tot_before, counts_after=tot_after,
	counts_removed_frac=0, stringsAsFactors=FALSE
	)

utils::write.table(tab, file=outcontam, quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)

cat(sprintf('[ambient] %s: method=%s, counts %.0f -> %.0f (%.2f%% removed), %d genes x %d cells\n',
	sampleid, method, tot_before, tot_after,
	100*(tot_before-tot_after)/tot_before, nrow(corrected), ncol(corrected)))
