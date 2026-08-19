# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:
"""Filter cells on total UMI, detected genes and mitochondrial percentage.

scanpy replacement for the Seurat implementation in v0.1.0. The three criteria
are unchanged, so thresholds carry over and cell counts are directly comparable:

	Seurat nCount_RNA        == scanpy total_counts
	Seurat nFeature_RNA      == scanpy n_genes_by_counts
	PercentageFeatureSet('^MT-|^mt-') == pct_counts_mt   (the default `mt` gene set)

Three criteria, but more than three metrics. Every gene set in `config.geneset`
gets a `pct_counts_<set>` column and a violin panel -- mitochondrial, ribosomal
and hemoglobin by default. Only the mitochondrial one has a threshold. The other
sets are recorded and plotted because they are how a reader recognises what kind
of bad a cell is (ribosomal-dominated dying cells, hemoglobin-dominated red
blood cell contamination), and left unfiltered because a defensible cut-off for
either is tissue-dependent -- cellqc does not exclude cells on a number nobody
has looked at yet.

Every one of them is computed on the *ambient-corrected* matrix, which is the
matrix every later stage and the user work with. The same metrics from the
uncorrected Cell Ranger counts are carried alongside under a `raw_` prefix, so
what the ambient correction took off a cell is visible in `.obs` -- and so the
filter thresholds can be reasoned about on both scales.

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
genesets = dict(snakemake.params['geneset'])

qcutil.set_seed(snakemake.params['seed'])

MITO_SET = qcutil.MITO_SET
# `mt` first, then the rest in config order: the mitochondrial panel is the one
# with a threshold and the one a reader looks at first.
SET_NAMES = [MITO_SET] + [n for n in genesets if n != MITO_SET]
PCT = {name: f'pct_counts_{name}' for name in SET_NAMES}

# The metrics recomputed on the uncorrected counts, and the prefix they get.
# `raw_` means *before ambient correction* -- the source is the Cell Ranger
# filtered matrix, the same cells as the corrected one, NOT the all-droplets
# raw_feature_bc_matrix.h5 that `barcoderank` reads.
RAW_METRICS = ('total_counts', 'n_genes_by_counts') + tuple(PCT[n] for n in SET_NAMES)
RAW_PREFIX = 'raw_'


def qc_metrics(adata):
	"""Flag every gene set in `.var` and let scanpy compute its percentages.

	Returns `(adata, matched)` where `matched` maps a set to `(n_genes,
	matched_by)`. The membership flags stay in `.var` -- which genes were counted
	is not recoverable from the percentage afterwards, and a reference that names
	its mitochondrial genes unusually is exactly the one whose numbers get
	questioned.
	"""
	matched = {}
	for name in SET_NAMES:
		mask, matched_by = qcutil.gene_set_mask(adata.var_names, genesets[name])
		adata.var[name] = mask
		matched[name] = (int(mask.sum()), matched_by)
	sc.pp.calculate_qc_metrics(
		adata, qc_vars=SET_NAMES, percent_top=None, log1p=False, inplace=True)
	return adata, matched


def add_raw_metrics(adata):
	"""Attach the pre-correction metrics as `raw_total_counts` & co.

	Without them the only per-cell numbers in the final `.obs` are post-ambient,
	and how much a given cell lost to the correction can be recovered only by
	reopening the Cell Ranger matrix and joining it by hand. Same metric
	definitions, same gene sets, so the pairs are comparable.
	"""
	raw = sc.read_10x_h5(rawfile)
	raw.var_names_make_unique()
	raw, _ = qc_metrics(raw)

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

	adata, matched = qc_metrics(adata)
	adata = add_raw_metrics(adata)
	n_before = adata.n_obs
	n_mt_genes = matched[MITO_SET][0]
	print(
		f'[filterbycount] {sampleid}: {n_before} cells x {adata.n_vars} genes',
		flush=True,
		)
	for name in SET_NAMES:
		n_genes, matched_by = matched[name]
		print(
			f'[filterbycount] {sampleid}: gene set {name!r}: {n_genes} genes '
			f'(matched by {matched_by}), median {np.nanmedian(adata.obs[PCT[name]]):.2f}%'
			+ ('' if name == MITO_SET else ' -- recorded, not filtered on'),
			flush=True,
			)
	if n_mt_genes == 0:
		# Only the mitochondrial set is worth an alarm: it is the one with a
		# threshold, so an empty set there means the filter silently does nothing.
		print(
			f'[filterbycount] {sampleid}: WARNING no gene matched the {MITO_SET!r} set '
			f'({genesets[MITO_SET]}). pct_counts_mt is 0 for every cell, so the mito '
			'threshold excludes nothing. Check that var_names are gene symbols for the '
			'expected organism, and add the reference\'s names to geneset.mt if not.',
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
		}])
	# Three columns per gene set: how many genes it matched, how they were
	# recognised, and the median percentage over the cells entering the filter
	# (pre-filter, like ncell_before). All three belong in the report -- a median
	# % is not comparable across cohorts without knowing which genes it was
	# computed over, or whether the reference matched by pattern or by fallback.
	for name in SET_NAMES:
		n_genes, matched_by = matched[name]
		stat[f'n_{name}_genes'] = n_genes
		stat[f'{name}_matched_by'] = matched_by
		stat[f'median_{PCT[name]}'] = round(float(np.nanmedian(adata.obs[PCT[name]])), 4)
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


def violin_panels():
	"""`(column, label, threshold, kind)` per panel, thresholds first.

	The two count criteria, then one panel per gene set. Only the mitochondrial
	set has a threshold; the others pass `None`, which is drawn as a panel with no
	dashed line and an explicit '(not filtered)' label rather than being left to
	look like a criterion whose line happens to be off-scale.
	"""
	panels = [
		('total_counts', 'Total UMI', mincount, 'min'),
		('n_genes_by_counts', 'Detected genes', minfeature, 'min'),
		]
	for name in SET_NAMES:
		label = qcutil.gene_set_label(name, genesets[name])
		if name == MITO_SET:
			panels.append((PCT[name], label, mito, 'max'))
		else:
			panels.append((PCT[name], f'{label}\n(not filtered)', None, None))
	return panels


def violin(adata, outfile, subtitle):
	"""QC violins with the applied thresholds drawn as lines."""
	qcutil.setup_matplotlib()
	import matplotlib.pyplot as plt

	panels = violin_panels()
	# At most three panels per row. One long row would keep the v0.3.2 look, but
	# the slide deck scales every figure to \linewidth with keepaspectratio, so a
	# 5-wide strip lands on the slide at ~60% of the previous size and its axis
	# labels stop being readable. Wrapping keeps the panels the size they were
	# (2.5 x 3.2in each, the v0.3.2 geometry) whatever a cohort's gene sets are;
	# three or fewer panels give exactly the v0.3.2 single-row figure.
	ncol = min(len(panels), 3)
	nrow = -(-len(panels) // ncol)
	fig, axes = plt.subplots(nrow, ncol, figsize=(2.5 * ncol, 3.2 * nrow), squeeze=False)
	axes = axes.ravel()
	for ax in axes[len(panels):]:
		ax.set_axis_off()
	for ax, (col, label, thr, kind) in zip(axes, panels):
		vals = adata.obs[col].to_numpy().astype(float)
		parts = ax.violinplot([vals], showextrema=False, widths=0.8)
		for body in parts['bodies']:
			body.set_facecolor('#106e78')
			body.set_alpha(0.55)
		ax.boxplot([vals], widths=0.12, showfliers=False,
			medianprops=dict(color='black'), whiskerprops=dict(color='black'),
			capprops=dict(color='black'), boxprops=dict(color='black'))
		if thr is not None:
			ax.axhline(thr, color='#c44e52', linestyle='--', linewidth=1)
			ax.annotate(f'{kind} {thr:g}', xy=(1.02, thr), xycoords=('axes fraction', 'data'),
				color='#c44e52', fontsize=7, va='center')
		ax.set_title(label, fontsize=9)
		ax.set_xticks([])
		if not col.startswith('pct_counts_'):
			ax.set_yscale('log')
	fig.suptitle(f'{sampleid} — {subtitle} (n={adata.n_obs})', fontsize=10)
	if nrow == 1:
		fig.tight_layout()
	else:
		# tight_layout does not reserve room for suptitle across stacked rows.
		fig.tight_layout(rect=(0, 0, 1, 0.98))
	qcutil.savefig(fig, outfile[:-len('.pdf')])
	plt.close(fig)


if __name__ == '__main__':
	main()
