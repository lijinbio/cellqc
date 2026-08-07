# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:
"""Filter cells on total UMI, detected genes and mitochondrial percentage.

scanpy replacement for the Seurat implementation in v0.1.0. The three criteria
and the mitochondrial gene pattern are unchanged, so thresholds carry over and
cell counts are directly comparable:

	Seurat nCount_RNA        == scanpy total_counts
	Seurat nFeature_RNA      == scanpy n_genes_by_counts
	PercentageFeatureSet('^MT-|^mt-') == pct_counts_mt

Those three are computed on the *ambient-corrected* matrix, which is the matrix
every later stage and the user work with. The same three from the uncorrected
Cell Ranger counts are carried alongside under a `raw_` prefix, so what the
ambient correction took off a cell is visible in `.obs` -- and so the filter
thresholds can be reasoned about on both scales.

v0.1.0 reported only the cell count before and after, which hides which
threshold did the work. Every criterion is counted separately here, including
overlaps, so an exclusion can always be attributed.
"""

import numpy as np
import pandas as pd
import scanpy as sc

from cellqc import qcutil

infile = snakemake.input['corrected']
rawfile = snakemake.input['raw']
out_h5ad = snakemake.output[0]
out_violin_before = snakemake.output[1]
out_violin_after = snakemake.output[3]
out_stat = snakemake.output[5]

mincount = snakemake.params['mincount']
minfeature = snakemake.params['minfeature']
mito = snakemake.params['mito']
sampleid = snakemake.params['sampleid']

qcutil.set_seed(snakemake.params['seed'])

MITO_PREFIXES = qcutil.MITO_PREFIXES

# The metrics recomputed on the uncorrected counts, and the prefix they get.
# `raw_` means *before ambient correction* -- the source is the Cell Ranger
# filtered matrix, the same cells as the corrected one, NOT the all-droplets
# raw_feature_bc_matrix.h5 that `barcoderank` reads.
RAW_METRICS = ('total_counts', 'n_genes_by_counts', 'pct_counts_mt')
RAW_PREFIX = 'raw_'


def qc_metrics(adata):
	adata.var['mt'] = [n.startswith(MITO_PREFIXES) for n in adata.var_names]
	sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
	return adata


def add_raw_metrics(adata):
	"""Attach the pre-correction metrics as `raw_total_counts` & co.

	Without them the only per-cell numbers in the final `.obs` are post-ambient,
	and how much a given cell lost to the correction can be recovered only by
	reopening the Cell Ranger matrix and joining it by hand. Same metric
	definitions, same mitochondrial pattern, so the pairs are comparable.
	"""
	raw = sc.read_10x_h5(rawfile)
	raw.var_names_make_unique()
	raw = qc_metrics(raw)

	missing = adata.obs.index.difference(raw.obs.index)
	if len(missing):
		raise ValueError(
			f'{sampleid}: {len(missing)} of {adata.n_obs} barcodes are absent from {rawfile} '
			f'(e.g. {list(missing[:3])}). Refusing to write {RAW_PREFIX}* columns that would be '
			'silently NaN.'
			)
	raw_obs = raw.obs.reindex(adata.obs.index)
	for col in RAW_METRICS:
		adata.obs[RAW_PREFIX + col] = raw_obs[col].to_numpy()
	return adata


def main():
	adata = sc.read_10x_h5(infile)
	adata.var_names_make_unique()
	adata.obs['sampleid'] = sampleid

	# Inspect before analysing: assert the matrix is what the rest of the
	# pipeline assumes (integer UMI counts), rather than discovering it later.
	x = adata.X
	sample = x.data[:1000] if hasattr(x, 'data') else np.asarray(x).ravel()[:1000]
	if sample.size and not np.allclose(sample, np.round(sample)):
		raise ValueError(
			f'{sampleid}: input matrix is not integer-valued. Every downstream '
			'count-based model (DoubletFinder, scDblFinder) assumes UMI counts.'
			)

	adata = qc_metrics(adata)
	adata = add_raw_metrics(adata)
	n_before = adata.n_obs
	n_mt_genes = int(adata.var['mt'].sum())
	print(
		f'[filterbycount] {sampleid}: {n_before} cells x {adata.n_vars} genes, '
		f'{n_mt_genes} mitochondrial genes matched by {MITO_PREFIXES}',
		flush=True,
		)
	if n_mt_genes == 0:
		print(
			f'[filterbycount] {sampleid}: WARNING no gene matched {MITO_PREFIXES}. '
			'pct_counts_mt is 0 for every cell, so the mito threshold excludes nothing. '
			'Check that var_names are gene symbols for the expected organism.',
			flush=True,
			)

	# What the corrected/uncorrected pair is for, stated once per sample: a
	# median well away from the reported contamination is worth looking at.
	raw_tot = adata.obs[RAW_PREFIX + 'total_counts'].to_numpy(dtype=float)
	tot = adata.obs['total_counts'].to_numpy(dtype=float)
	removed = np.divide(raw_tot - tot, raw_tot, out=np.full(raw_tot.shape, np.nan), where=raw_tot > 0)
	print(
		f'[filterbycount] {sampleid}: ambient correction removed a median '
		f'{100 * np.nanmedian(removed):.2f}% of a cell\'s UMI '
		f'({RAW_PREFIX}* columns carry the pre-correction metrics)',
		flush=True,
		)

	violin(adata, out_violin_before, 'before filtering')

	fail_count = adata.obs['total_counts'].to_numpy() < mincount
	fail_feature = adata.obs['n_genes_by_counts'].to_numpy() < minfeature
	fail_mito = adata.obs['pct_counts_mt'].to_numpy() > mito
	keep = ~(fail_count | fail_feature | fail_mito)

	# Attribute every exclusion, including cells failing more than one criterion.
	stat = pd.DataFrame([{
		'sampleid': sampleid,
		'ncell_before': int(n_before),
		'ncell_after': int(keep.sum()),
		'ncell_removed': int((~keep).sum()),
		'frac_removed': float((~keep).sum() / n_before) if n_before else np.nan,
		'fail_mincount': int(fail_count.sum()),
		'fail_minfeature': int(fail_feature.sum()),
		'fail_mito': int(fail_mito.sum()),
		'fail_mincount_only': int((fail_count & ~fail_feature & ~fail_mito).sum()),
		'fail_minfeature_only': int((fail_feature & ~fail_count & ~fail_mito).sum()),
		'fail_mito_only': int((fail_mito & ~fail_count & ~fail_feature).sum()),
		'fail_multiple': int(((fail_count.astype(int) + fail_feature.astype(int) + fail_mito.astype(int)) > 1).sum()),
		'mincount': mincount,
		'minfeature': minfeature,
		'mito': mito,
		'n_mito_genes': n_mt_genes,
		}])
	stat.to_csv(out_stat, sep='\t', index=False)

	print(
		f'[filterbycount] {sampleid}: {n_before} -> {int(keep.sum())} cells '
		f'({int((~keep).sum())} removed: {int(fail_count.sum())} by mincount>={mincount}, '
		f'{int(fail_feature.sum())} by minfeature>={minfeature}, '
		f'{int(fail_mito.sum())} by mito<={mito}%; criteria overlap)',
		flush=True,
		)
	if keep.sum() == 0:
		raise ValueError(
			f'{sampleid}: all {n_before} cells were removed by filterbycount. '
			f'Thresholds (mincount={mincount}, minfeature={minfeature}, mito={mito}) '
			'are almost certainly wrong for this sample.'
			)

	adata = adata[keep].copy()
	violin(adata, out_violin_after, 'after filtering')
	adata.write_h5ad(out_h5ad, compression='gzip')


def violin(adata, outfile, subtitle):
	"""QC violins with the applied thresholds drawn as lines."""
	qcutil.setup_matplotlib()
	import matplotlib.pyplot as plt

	panels = [
		('total_counts', 'Total UMI', mincount, 'min'),
		('n_genes_by_counts', 'Detected genes', minfeature, 'min'),
		('pct_counts_mt', '% mitochondrial', mito, 'max'),
		]
	fig, axes = plt.subplots(1, 3, figsize=(7.5, 3.2))
	for ax, (col, label, thr, kind) in zip(axes, panels):
		vals = adata.obs[col].to_numpy().astype(float)
		parts = ax.violinplot([vals], showextrema=False, widths=0.8)
		for body in parts['bodies']:
			body.set_facecolor('#106e78')
			body.set_alpha(0.55)
		ax.boxplot([vals], widths=0.12, showfliers=False,
			medianprops=dict(color='black'), whiskerprops=dict(color='black'),
			capprops=dict(color='black'), boxprops=dict(color='black'))
		ax.axhline(thr, color='#c44e52', linestyle='--', linewidth=1)
		ax.annotate(f'{kind} {thr:g}', xy=(1.02, thr), xycoords=('axes fraction', 'data'),
			color='#c44e52', fontsize=7, va='center')
		ax.set_title(label, fontsize=9)
		ax.set_xticks([])
		if col != 'pct_counts_mt':
			ax.set_yscale('log')
	fig.suptitle(f'{sampleid} — {subtitle} (n={adata.n_obs})', fontsize=10)
	fig.tight_layout()
	qcutil.savefig(fig, outfile[:-len('.pdf')])
	plt.close(fig)


if __name__ == '__main__':
	main()
