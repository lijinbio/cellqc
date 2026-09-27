# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:
"""Apply the doublet calls and write the QC'd matrix.

Every caller listed in `doublet.run` contributes its score and class to .obs
under a namespaced column, but only `doublet.decider` removes cells. Keeping the
decision with one caller avoids an undeclared ensemble rule: a union of two
callers removes more cells than the assumed multiplet rate, an intersection
fewer, and either is a choice that should be made deliberately rather than
emerge from how the code happens to be wired.

Concordance between callers is computed here (2x2 table plus Cohen's kappa) and
carried into both reports. It is a consistency measure, not an accuracy measure
-- with no ground-truth doublets neither caller can be shown correct.

A caller that failed on this sample is tolerated unless it is the decider: its
columns are absent from .obs, its summary row is NA, and no concordance is
computed against it. A failed decider skips this step -- and so excludes the
sample -- rather than silently handing the decision to another caller, which
would make one sample's doublet removal differ in method from the cohort's.
"""

import itertools

import numpy as np
import pandas as pd
import scanpy as sc

from cellqc import qcstatus, qcutil

in_h5ad = snakemake.input['h5ad']
in_meta = qcstatus.as_list(snakemake.input['metadata'])
in_caller_status = qcstatus.as_list(snakemake.input['caller_status'])
out_h5ad = snakemake.output['h5ad']
out_summary = snakemake.output['summary']
out_concordance = snakemake.output['concordance']

CONCORDANCE_COLUMNS = [
	'sampleid', 'caller_a', 'caller_b', 'both_doublet', 'only_a', 'only_b',
	'both_singlet', 'kappa',
	]

sampleid = snakemake.params['sampleid']
callers = list(snakemake.params['callers'])
decider = snakemake.params['decider']

# There is no skip flag: the schema requires at least one caller in doublet.run
# and config.smk requires the decider to be one of them. An empty list here
# would mean the DAG was built from a config the validator should have rejected.
if not callers:
	raise ValueError(f'{sampleid}: no doublet caller configured; doublet.run must list at least one.')


def cohens_kappa(a, b):
	"""Chance-corrected agreement between two binary label vectors."""
	a = np.asarray(a)
	b = np.asarray(b)
	n = len(a)
	if n == 0:
		return np.nan
	po = float((a == b).mean())
	pa1, pb1 = float(a.mean()), float(b.mean())
	pe = pa1 * pb1 + (1 - pa1) * (1 - pb1)
	if np.isclose(pe, 1.0):
		return np.nan
	return (po - pe) / (1 - pe)


def main(st):
	adata = sc.read_h5ad(in_h5ad)
	n_before = adata.n_obs
	rows = []

	# `metadata` and `caller_status` are both in `callers` order (common.smk).
	ran = [c for c, sf in zip(callers, in_caller_status) if qcstatus.usable(sf)]
	lost = [c for c in callers if c not in ran]
	if lost:
		st.note(f"{', '.join(lost)} unavailable; {decider} alone decided and no concordance was computed")

	# Merge each caller's metadata onto .obs by barcode.
	for caller, path in zip(callers, in_meta):
		if caller not in ran:
			continue
		meta = pd.read_csv(path, sep='\t', header=0).set_index('barcode')
		missing = adata.obs_names.difference(meta.index)
		if len(missing):
			raise ValueError(
				f'{sampleid}: {caller} returned no call for {len(missing)} of {adata.n_obs} '
				'barcodes. Refusing to write a matrix with silently missing doublet calls.'
				)
		meta = meta.reindex(adata.obs_names)
		for col in meta.columns:
			adata.obs[col] = meta[col].to_numpy()

	for caller in callers:
		if caller not in ran:
			rows.append({
				'sampleid': sampleid, 'caller': caller, 'is_decider': False,
				'ndoublet': np.nan, 'ncell': n_before, 'frac_doublet': np.nan,
				})
			print(f'[filterdoublet] {sampleid}: {caller} did not run on this sample', flush=True)
			continue
		_, classcol = qcutil.CALLER_COLUMNS[caller]
		if classcol not in adata.obs:
			raise ValueError(f'{sampleid}: expected column {classcol} from caller {caller}, not found')
		adata.obs[classcol] = pd.Categorical(adata.obs[classcol], categories=['Singlet', 'Doublet'])
		nd = int((adata.obs[classcol] == 'Doublet').sum())
		rows.append({
			'sampleid': sampleid, 'caller': caller, 'is_decider': caller == decider,
			'ndoublet': nd, 'ncell': n_before, 'frac_doublet': nd / n_before if n_before else np.nan,
			})
		print(f'[filterdoublet] {sampleid}: {caller} called {nd}/{n_before} doublets '
			f'({100 * nd / n_before:.2f}%)' + ('  <- decides removal' if caller == decider else '  (diagnostic only)'),
			flush=True)

	summary = pd.DataFrame(rows)

	# Pairwise concordance. Descriptive: agreement is not evidence of accuracy.
	conc = []
	for a, b in itertools.combinations(ran, 2):
		ca = (adata.obs[qcutil.CALLER_COLUMNS[a][1]] == 'Doublet').to_numpy()
		cb = (adata.obs[qcutil.CALLER_COLUMNS[b][1]] == 'Doublet').to_numpy()
		k = cohens_kappa(ca, cb)
		conc.append({
			'sampleid': sampleid, 'caller_a': a, 'caller_b': b,
			'both_doublet': int((ca & cb).sum()),
			'only_a': int((ca & ~cb).sum()),
			'only_b': int((~ca & cb).sum()),
			'both_singlet': int((~ca & ~cb).sum()),
			'kappa': k,
			})
		note = ''
		if np.isfinite(k) and k < 0.5:
			note = ('  WARNING low agreement: the doublet call is unstable for this '
				'sample and should not be treated as settled')
		print(f'[filterdoublet] {sampleid}: {a} vs {b} Cohen kappa={k:.3f} '
			f'(both={int((ca & cb).sum())}, only_{a}={int((ca & ~cb).sum())}, '
			f'only_{b}={int((~ca & cb).sum())}){note}', flush=True)

	_, decidercol = qcutil.CALLER_COLUMNS[decider]
	keep = (adata.obs[decidercol] != 'Doublet').to_numpy()
	ndoublet = int((~keep).sum())

	if keep.sum() == 0:
		raise ValueError(f'{sampleid}: {decider} classified all {n_before} cells as doublets.')

	adata.uns['cellqc_doublet_decider'] = decider
	adata.uns['cellqc_doublet_callers'] = list(callers)
	adata.uns['cellqc_homotypic_modelled'] = False

	filtered = adata[keep].copy()
	filtered.write_h5ad(out_h5ad, compression='gzip')

	summary['ncell_before'] = n_before
	summary['ncell_after'] = int(keep.sum())
	summary['decider'] = decider
	summary.to_csv(out_summary, sep='\t', index=False)

	# Always written, even with a single caller, so the reports can read it
	# unconditionally.
	cdf = pd.DataFrame(conc, columns=CONCORDANCE_COLUMNS) if conc else pd.DataFrame(columns=CONCORDANCE_COLUMNS)
	cdf.to_csv(out_concordance, sep='\t', index=False)

	print(f'[filterdoublet] {sampleid}: {n_before} -> {int(keep.sum())} cells '
		f'({ndoublet} removed by {decider}; homotypic doublets not modelled, so this '
		'is a slight over-removal)', flush=True)


if __name__ == '__main__':
	qcstatus.run(snakemake, 'filterdoublet', main, requires=snakemake.input['upstream'])
