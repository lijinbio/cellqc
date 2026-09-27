# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:
"""Collect every sample's step statuses into the cohort's QC status.

Writes

	result/qc_status.csv  sample, step, status, message -- every sample, every step
	result/manifest.tsv   sample, included, ncell, reason, fallback -- one row each

A sample is included in result/ when its last step, postproc, is usable. That
single test is enough because each required step skips when the one before it
is unusable (ambient -> filterbycount -> doublet decider -> filterdoublet ->
postproc), and the skip message carries the root cause forward, so `reason`
names the step that actually went wrong. Diagnostic steps -- barcoderank, the
nuclear fraction, a non-deciding doublet caller, a comparison-only ambient
method -- are reported but never exclude a sample.

Every sample gets a final `result` row, `ok` or `skipped`, so the status table
alone says which samples reached result/.
"""

import os
import sys

import pandas as pd

from cellqc import qcstatus

status_files = qcstatus.as_list(snakemake.input['status'])
summary_files = qcstatus.as_list(snakemake.input['summary'])
ids = list(snakemake.params['samples'])
steps = dict(snakemake.params['steps'])
has_bam = dict(snakemake.params['has_bam'])
out_csv = snakemake.output['csv']
out_manifest = snakemake.output['manifest']

RESULT_STEP = 'result'


def status_path(step, sample):
	return f'{step}/{sample}_status.tsv'


def ncell_final(sample):
	path = f'filterdoublet/{sample}_doublet_summary.txt'
	if not os.path.exists(path) or os.path.getsize(path) == 0:
		return None
	df = pd.read_csv(path, sep='\t')
	return int(df['ncell_after'].iloc[0]) if 'ncell_after' in df and len(df) else None


def main():
	known = set(status_files)
	rows, manifest = [], []
	for s in ids:
		per = []
		for step in steps[s]:
			path = status_path(step, s)
			got = qcstatus.read_status(path) if path in known else None
			if got is None:
				got = [{'sample': s, 'step': step, 'status': qcstatus.FAILED, 'message': f'no status recorded in {path}'}]
			per += got
			if step == 'ambient' and not has_bam[s]:
				# Not a rule that ran and skipped -- it is never scheduled without a
				# BAM -- but a reader of this table should not have to know that.
				per.append({'sample': s, 'step': 'nuclear_fraction', 'status': qcstatus.SKIPPED,
					'message': 'no indexed possorted_genome_bam.bam in the Cell Ranger directory'})

		last = next(r for r in per if r['step'] == 'postproc')
		included = last['status'] in qcstatus.USABLE
		ncell = ncell_final(s) if included else None
		reason = '' if included else qcstatus.reason(last)
		per.append({
			'sample': s, 'step': RESULT_STEP,
			'status': qcstatus.OK if included else qcstatus.SKIPPED,
			'message': f'included ({ncell} cells)' if included else f'excluded: {reason}',
			})
		fallbacks = '; '.join(f"{r['step']}: {r['message']}" for r in per if r['status'] == qcstatus.FALLBACK)
		manifest.append({'sample': s, 'included': included, 'ncell': ncell, 'reason': reason, 'fallback': fallbacks})
		rows += per

	table = pd.DataFrame(rows, columns=list(qcstatus.COLUMNS))
	table.to_csv(out_csv, index=False)
	mf = pd.DataFrame(manifest, columns=['sample', 'included', 'ncell', 'reason', 'fallback'])
	mf['ncell'] = mf['ncell'].astype('Int64')
	mf.to_csv(out_manifest, sep='\t', index=False, na_rep='NA')

	n_inc = int(mf['included'].sum())
	counts = table[table['step'] != RESULT_STEP]['status'].value_counts().to_dict()
	print(f'[qcstatus] {n_inc}/{len(mf)} sample(s) included in result/; step outcomes: {counts}', flush=True)
	for _, r in table[~table['status'].isin([qcstatus.OK])].iterrows():
		if r['step'] != RESULT_STEP:
			print(f"[qcstatus] {r['sample']}: {r['step']} {r['status']}: {r['message']}", file=sys.stderr, flush=True)
	for _, r in mf[~mf['included']].iterrows():
		print(f"[qcstatus] EXCLUDED {r['sample']}: {r['reason']}", file=sys.stderr, flush=True)


if __name__ == '__main__':
	main()
