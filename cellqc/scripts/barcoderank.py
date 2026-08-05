# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:
"""Barcode rank ("knee") plot from the raw Cell Ranger matrix.

Purely diagnostic: cellqc does not call cells, it uses Cell Ranger's EmptyDrops
call. The plot exists so the separation between cells and empty droplets can be
judged by eye, and so a library where that separation is poor is visible rather
than inferred from a downstream oddity.

The knee and inflection are computed the way DropletUtils::barcodeRanks does --
on a spline fitted to the log-log curve over the region between the total-count
range -- and are reported as descriptive annotations only. Nothing downstream
uses them.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.interpolate import LSQUnivariateSpline

from cellqc import qcutil

in_raw = snakemake.input['raw']
in_filtered = snakemake.input['filtered']
out_pdf, out_png, out_knee = snakemake.output[0], snakemake.output[1], snakemake.output[2]
sampleid = snakemake.params['sampleid']


def totals_from_h5(path):
	adata = sc.read_10x_h5(path)
	return np.asarray(adata.X.sum(axis=1)).ravel(), adata.n_obs


def knee_inflection(totals, lower=100, df=20):
	"""Knee and inflection on the log-log rank/total curve.

	Follows DropletUtils::barcodeRanks: run-length encode the sorted totals so
	each DISTINCT total count contributes one point at its mean rank, keep those
	above `lower`, fit a spline to log(total) vs log(rank), then take

	  inflection = steepest descent  (minimum first derivative)
	  knee       = sharpest bend     (minimum signed curvature)

	Two details matter and both were got wrong on the first attempt, so they are
	recorded here:

	  * Ranking per barcode rather than per distinct total lets the huge ambient
	    plateau (~1.5e6 barcodes at 1-100 UMI) dominate the fit, which pushes the
	    derivative minimum onto the `lower` cliff and yields a knee of ~100.
	  * Evenly spaced knots do the same thing, because point density along x is
	    wildly uneven. Knots go at quantiles of x.

	Returns (knee_total, inflection_total, n_above_lower); NaN when the curve is
	too short to fit, reported rather than guessed.
	"""
	o = np.sort(totals)[::-1]
	vals, counts = np.unique(o[::-1], return_counts=True)
	vals, counts = vals[::-1], counts[::-1]
	# Mean rank within each run of equal totals, as barcodeRanks does.
	run_rank = np.cumsum(counts) - (counts - 1) / 2

	keep = vals > lower
	n = int(counts[keep].sum())
	if keep.sum() < 50:
		return np.nan, np.nan, n

	x = np.log10(run_rank[keep])
	y = np.log10(vals[keep].astype(float))
	knots = np.unique(np.quantile(x, np.linspace(0, 1, df + 2)[1:-1]))
	if len(knots) < 3:
		return np.nan, np.nan, n
	try:
		spl = LSQUnivariateSpline(x, y, knots, k=3)
	except Exception:
		return np.nan, np.nan, n

	d1 = spl.derivative(1)(x)
	d2 = spl.derivative(2)(x)
	curvature = d2 / np.power(1 + d1 ** 2, 1.5)
	i_inflect = int(np.argmin(d1))
	i_knee = int(np.argmin(curvature))
	return float(10 ** y[i_knee]), float(10 ** y[i_inflect]), n


def main():
	raw_totals, n_raw = totals_from_h5(in_raw)
	_, n_called = totals_from_h5(in_filtered)

	nonzero = raw_totals[raw_totals > 0]
	print(f'[barcoderank] {sampleid}: {n_raw} raw barcodes, {len(nonzero)} with >0 UMI, '
		f'{n_called} called as cells by Cell Ranger', flush=True)

	knee, inflection, n_above = knee_inflection(raw_totals)
	order = np.sort(nonzero)[::-1]
	rank = np.arange(1, len(order) + 1)
	# UMI threshold implied by Cell Ranger's call, i.e. the total of the lowest
	# ranked called cell -- useful next to the knee for judging the call.
	cr_threshold = float(np.sort(raw_totals)[::-1][n_called - 1]) if n_called <= len(raw_totals) else np.nan

	pd.DataFrame([{
		'sampleid': sampleid,
		'n_raw_barcodes': int(n_raw),
		'n_nonzero_barcodes': int(len(nonzero)),
		'n_called_cells': int(n_called),
		'cellranger_umi_threshold': cr_threshold,
		'knee_total': knee,
		'inflection_total': inflection,
		'n_barcodes_above_100umi': int(n_above),
		}]).to_csv(out_knee, sep='\t', index=False)

	qcutil.setup_matplotlib()
	import matplotlib.pyplot as plt

	fig, ax = plt.subplots(figsize=(4.6, 4.0))
	# Rasterized inside the vector PDF: ~1e6 points would otherwise make a file
	# no viewer can open, while the axes and labels stay editable text.
	ax.plot(rank, order, color='#106e78', linewidth=1.2, rasterized=True)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.axvline(n_called, color='#c44e52', linestyle='--', linewidth=1)
	ax.annotate(f'Cell Ranger: {n_called:,} cells', xy=(n_called, ax.get_ylim()[1]),
		xytext=(4, -4), textcoords='offset points', color='#c44e52', fontsize=7,
		ha='left', va='top', rotation=90)
	for val, label, color in ((knee, 'knee', '#4c72b0'), (inflection, 'inflection', '#8172b2')):
		if np.isfinite(val):
			ax.axhline(val, color=color, linestyle=':', linewidth=1)
			ax.annotate(f'{label} {val:,.0f}', xy=(1.0, val), xycoords=('axes fraction', 'data'),
				xytext=(-2, 2), textcoords='offset points', color=color, fontsize=7, ha='right')
	ax.set_xlabel('Barcode rank')
	ax.set_ylabel('Total UMI')
	ax.set_title(f'{sampleid}', fontsize=10)
	fig.tight_layout()
	qcutil.savefig(fig, out_pdf[:-len('.pdf')])
	plt.close(fig)

	print(f'[barcoderank] {sampleid}: knee={knee}, inflection={inflection}, '
		f'Cell Ranger UMI threshold={cr_threshold}', flush=True)


if __name__ == '__main__':
	main()
