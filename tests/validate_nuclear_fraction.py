#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:
"""Acceptance gate for the pysam nuclear-fraction reimplementation.

v0.2.0 replaces DropletQC::nuclear_fraction_tags() with a pysam single pass, to
remove a GitHub-only R dependency. That is only acceptable if it reproduces the
reference. This compares the two on the same sample and FAILS if it does not.

Gate (docs/design.md 2.2):
	identical barcode set
	Pearson r and Spearman rho > 0.999
	median |delta| < 0.001 and max |delta| < 0.01

If the gate fails the correct action is to revert to DropletQC, not to loosen
the threshold: matching the reference matters more than removing one install.

Usage:
	validate_nuclear_fraction.py <new.txt.gz> <dropletqc_reference.txt.gz> [outdir]
"""

import sys

import numpy as np
import pandas as pd
from scipy import stats

PEARSON_MIN = 0.999
SPEARMAN_MIN = 0.999
MEDIAN_ABS_MAX = 0.001
MAX_ABS_MAX = 0.01


def main(newfile, reffile, outdir=None):
	new = pd.read_csv(newfile, sep='\t').set_index('barcode')
	ref = pd.read_csv(reffile, sep='\t').set_index('barcode')

	print(f'new (pysam):     {len(new):,} barcodes, columns={list(new.columns)}')
	print(f'ref (DropletQC): {len(ref):,} barcodes, columns={list(ref.columns)}')

	failures = []

	only_new = new.index.difference(ref.index)
	only_ref = ref.index.difference(new.index)
	if len(only_new) or len(only_ref):
		failures.append(
			f'barcode sets differ: {len(only_new)} only in new, {len(only_ref)} only in reference')
	common = new.index.intersection(ref.index)
	print(f'shared barcodes: {len(common):,}')
	if len(common) == 0:
		print('FAIL: no shared barcodes'); return 1

	a = new.loc[common, 'nuclear_fraction'].to_numpy(dtype=float)
	b = ref.loc[common, 'nuclear_fraction'].to_numpy(dtype=float)
	ok = np.isfinite(a) & np.isfinite(b)
	n_nan = int((~ok).sum())
	if n_nan:
		print(f'NOTE: {n_nan} barcodes have NaN in one or both; excluded from correlation '
			f'(reported, not silently dropped)')
	a, b = a[ok], b[ok]

	pearson = stats.pearsonr(a, b)[0]
	spearman = stats.spearmanr(a, b)[0]
	delta = a - b
	med_abs = float(np.median(np.abs(delta)))
	max_abs = float(np.max(np.abs(delta)))

	print()
	print(f'  Pearson r        = {pearson:.6f}   (gate > {PEARSON_MIN})')
	print(f'  Spearman rho     = {spearman:.6f}   (gate > {SPEARMAN_MIN})')
	print(f'  median |delta|   = {med_abs:.6f}   (gate < {MEDIAN_ABS_MAX})')
	print(f'  max |delta|      = {max_abs:.6f}   (gate < {MAX_ABS_MAX})')
	print(f'  mean delta       = {float(np.mean(delta)):+.6f}  (systematic offset)')
	print(f'  new  median NF   = {float(np.median(a)):.4f}')
	print(f'  ref  median NF   = {float(np.median(b)):.4f}')

	if pearson <= PEARSON_MIN:
		failures.append(f'Pearson r={pearson:.6f} <= {PEARSON_MIN}')
	if spearman <= SPEARMAN_MIN:
		failures.append(f'Spearman rho={spearman:.6f} <= {SPEARMAN_MIN}')
	if med_abs >= MEDIAN_ABS_MAX:
		failures.append(f'median |delta|={med_abs:.6f} >= {MEDIAN_ABS_MAX}')
	if max_abs >= MAX_ABS_MAX:
		failures.append(f'max |delta|={max_abs:.6f} >= {MAX_ABS_MAX}')

	# Sensitivity of the statistic to counting secondary/duplicate alignments.
	if 'nuclear_fraction_primary_only' in new.columns:
		p = new.loc[common, 'nuclear_fraction_primary_only'].to_numpy(dtype=float)[ok]
		d = np.abs(a - p)
		print(f'\n  read-filter sensitivity: max |NF - NF_primary_only| = {np.nanmax(d):.6f}, '
			f'median = {np.nanmedian(d):.6f}')
		print('  (DropletQC counts secondary/duplicate alignments by default; matched here)')

	if outdir:
		plot(a, b, outdir)

	print()
	if failures:
		print('GATE FAILED:')
		for f in failures:
			print(f'  - {f}')
		print('\nAction: revert to DropletQC rather than loosening the gate.')
		return 1
	print('GATE PASSED: the pysam implementation reproduces DropletQC on this sample.')
	return 0


def plot(a, b, outdir):
	import os
	sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
	from cellqc import qcutil

	qcutil.setup_matplotlib()
	import matplotlib.pyplot as plt

	mean = (a + b) / 2
	diff = a - b
	fig, axes = plt.subplots(1, 2, figsize=(8, 3.6))
	axes[0].scatter(b, a, s=2, alpha=0.3, linewidths=0, color='#106e78', rasterized=True)
	lims = [min(a.min(), b.min()), max(a.max(), b.max())]
	axes[0].plot(lims, lims, color='#c44e52', linewidth=0.8, linestyle='--')
	axes[0].set_xlabel('DropletQC nuclear fraction')
	axes[0].set_ylabel('pysam nuclear fraction')
	axes[0].set_title('Agreement', fontsize=9)

	axes[1].scatter(mean, diff, s=2, alpha=0.3, linewidths=0, color='#106e78', rasterized=True)
	axes[1].axhline(0, color='#c44e52', linewidth=0.8, linestyle='--')
	axes[1].axhline(float(np.mean(diff)), color='#4c72b0', linewidth=0.8)
	axes[1].set_xlabel('mean of the two')
	axes[1].set_ylabel('pysam - DropletQC')
	axes[1].set_title('Bland-Altman', fontsize=9)
	fig.tight_layout()
	qcutil.savefig(fig, os.path.join(outdir, 'nf_validation'))
	print(f'wrote {outdir}/nf_validation.pdf')


if __name__ == '__main__':
	if len(sys.argv) < 3:
		print(__doc__)
		sys.exit(2)
	sys.exit(main(*sys.argv[1:4]))
