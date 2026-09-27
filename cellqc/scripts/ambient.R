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

source(snakemake@params[['guard']])

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
	png(outpng, width=6, height=4.5, units='in', res=300); fn(); dev.off()
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

# The figure for a sample whose applied method failed: the slot in both reports
# says what happened instead of showing nothing, or a stale figure.
plot_failed=function(m, why) {
	plot.new()
	title(main=sprintf('%s: %s failed', sampleid, m))
	text(0.5, 0.62, paste(strwrap(why, width=70), collapse='\n'), cex=0.75)
	text(0.5, 0.25, 'Uncorrected Cell Ranger counts were used downstream.', cex=0.85, font=2)
}

runners=list(soupx=run_soupx, decontx=run_decontx)

methods=unique(c(method, compare))
methods=methods[methods!='none']

cellqc_guard('ambient', {
	features=feature_table()

	# Each method is run on its own. A method that fails on this sample -- SoupX's
	# autoEstCont finds no marker genes in a low-complexity channel -- is recorded
	# and the others still run. If it is the applied method, the sample falls back
	# to the uncorrected counts (the `none` path) rather than being lost: ambient
	# correction improves a matrix, it is not what makes one usable.
	estimates=list()
	failures=list()
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
		est=tryCatch(runners[[m]](), error=function(e) e)
		if (inherits(est, 'error')) {
			failures[[m]]=cellqc_clean(conditionMessage(est))
			cat(sprintf('[ambient] %s: %s FAILED: %s\n', sampleid, m, failures[[m]]), file=stderr())
			if (m!=method) cellqc_substep(paste0('ambient.', m), 'failed', failures[[m]])
			next
		}
		estimates[[m]]=est
		if (m==method) corrected=est$counts
		else cellqc_substep(paste0('ambient.', m), 'ok', 'comparison only, not applied')
	}

	# The method whose counts were actually written: `method`, or `none` after a
	# fallback. The contamination table's `applied` column follows it.
	applied=method
	if (method=='none') {
		corrected=align_features(counts(read10xCounts(filt_h5, col.names=TRUE)), 'id')
		draw(function() { plot.new(); title(main=sprintf('%s: no ambient correction applied', sampleid)) })
	} else if (is.null(corrected)) {
		applied='none'
		corrected=align_features(counts(read10xCounts(filt_h5, col.names=TRUE)), 'id')
		draw(function() plot_failed(method, failures[[method]]))
		cellqc_fallback(sprintf('%s failed (%s); uncorrected Cell Ranger counts used', method, failures[[method]]))
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
	# side by side and can never be confused for one another. A method that failed
	# keeps its row, with NA estimates, so a failure is a visible gap rather than a
	# missing row; after a fallback a `none` row records what was applied instead.
	na_row=function(m, is_applied, ncell, after) data.frame(
		sampleid=sampleid, method=m, applied=is_applied,
		contamination_mean=NA_real_, contamination_median=NA_real_,
		contamination_min=NA_real_, contamination_max=NA_real_,
		ncell=ncell, counts_before=tot_before,
		counts_after=if (is_applied) after else NA_real_,
		counts_removed_frac=if (is_applied) (tot_before-after)/tot_before else NA_real_,
		stringsAsFactors=FALSE
		)
	rows=lapply(methods, function(m) {
		e=estimates[[m]]
		if (is.null(e)) return(na_row(m, FALSE, NA_integer_, NA_real_))
		data.frame(
			sampleid=sampleid, method=m, applied=(m==applied),
			contamination_mean=mean(e$contamination),
			contamination_median=stats::median(e$contamination),
			contamination_min=min(e$contamination),
			contamination_max=max(e$contamination),
			ncell=length(e$contamination),
			counts_before=tot_before,
			counts_after=ifelse(m==applied, tot_after, NA_real_),
			counts_removed_frac=ifelse(m==applied, (tot_before-tot_after)/tot_before, NA_real_),
			stringsAsFactors=FALSE
			)
		})
	# `method: none` with comparison methods keeps its v0.3.3 table (no `none` row);
	# with none at all, the `none` row is the whole table, as before.
	if (applied=='none' && (method!='none' || !length(rows)))
		rows=c(rows, list(na_row('none', TRUE, ncol(corrected), tot_after)))
	tab=do.call(rbind, rows)

	utils::write.table(tab, file=outcontam, quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)

	cat(sprintf('[ambient] %s: method=%s, counts %.0f -> %.0f (%.2f%% removed), %d genes x %d cells\n',
		sampleid, applied, tot_before, tot_after,
		100*(tot_before-tot_after)/tot_before, nrow(corrected), ncol(corrected)))
})
