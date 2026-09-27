# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:
"""Assemble everything both reports show, from the files the workflow wrote.

The HTML report and the PDF slide deck call this and nothing else, so the two
cannot drift into describing different runs. Nothing here recomputes an
analysis: every number is read back from the stage that produced it.

Every sample is in every table, including one that failed a step: what was
computed is shown, and what was not is NA rather than a missing row. A stats
file or figure that a failed or skipped step left as a 0-byte placeholder
(cellqc/qcstatus.py) is treated as not computed, and a figure that is missing
carries the reason from result/qc_status.csv instead. Nothing a missing input
does here may raise: the reports exist precisely to explain a run in which
something went wrong.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd

from cellqc import qcstatus

# The step whose status explains a missing figure.
FIGURE_STEP = {
	'barcoderank': 'barcoderank',
	'ambient': 'ambient',
	'violin_before': 'filterbycount',
	'violin_after': 'filterbycount',
	'nf': 'nuclear_fraction',
	'doubletfinder_pANN': 'doubletfinder',
	'doubletfinder_umap': 'doubletfinder',
	'scdblfinder_score': 'scdblfinder',
}


def computed(path):
	"""True for a file a step actually wrote; False for absent or a placeholder."""
	return os.path.exists(path) and os.path.getsize(path) > 0


def _read_concat(paths):
	frames = []
	for p in paths:
		if not computed(p):
			continue
		try:
			frames.append(pd.read_csv(p, sep='\t', header=0))
		except Exception as e:  # a damaged stats file must not take the report with it
			print(f'[reportdata] WARNING could not read {p}: {e}', flush=True)
	frames = [f for f in frames if len(f)]
	return pd.concat(frames, ignore_index=True) if frames else None


def _every_sample(table, ids):
	"""`table` with an all-NA row for each sample it lacks, in sample order.

	A sample that failed a step keeps its row in every table, so a gap reads as
	'not computed' instead of looking like a sample that was never in the run.
	"""
	if table is None or 'sampleid' not in table.columns:
		return table
	have = set(table['sampleid'].astype(str))
	gaps = [s for s in ids if s not in have]
	if gaps:
		# Keep integer and boolean columns integer and boolean: a NaN would turn
		# a column of cell counts into floats and print them as 1.022e+04.
		kinds = {c: table[c].dtype.kind for c in table.columns}
		table = pd.concat([table, pd.DataFrame({'sampleid': gaps})], ignore_index=True)
		for c, kind in kinds.items():
			if kind in 'iu':
				table[c] = table[c].astype('Int64')
			elif kind == 'b':
				table[c] = table[c].astype('boolean')
	order = {s: i for i, s in enumerate(ids)}
	key = table['sampleid'].map(order).fillna(len(ids))
	return table.iloc[np.argsort(key.to_numpy(), kind='stable')].reset_index(drop=True)


def read_status(path):
	"""result/qc_status.csv, or None if it is not there."""
	if not path or not computed(path):
		return None
	return pd.read_csv(path, keep_default_na=False, dtype=str)


def step_status(data, sid, step):
	"""The status row of one sample's step, as a dict, or None."""
	st = data.get('status')
	if st is None:
		return None
	m = st[(st['sample'] == sid) & (st['step'] == step)]
	return m.iloc[0].to_dict() if len(m) else None


def figure(data, key, sid, ext):
	"""`(path, None)` for a figure that exists, `(None, why)` for one that does not.

	`why` is the failing step's reason from the status table, so the report slot
	says e.g. `filterbycount failed: ...` rather than silently going missing.
	`(None, None)` means the figure was never expected for this sample (no BAM,
	no nuclear-fraction plot) and the slot is left out, as before.
	"""
	pattern = data['figures'].get(key, {}).get(sid)
	if pattern is None:
		return None, None
	path = pattern.format(ext=ext)
	if computed(path):
		return path, None
	row = step_status(data, sid, FIGURE_STEP.get(key, key))
	if row is not None and row['status'] not in qcstatus.USABLE:
		return None, qcstatus.reason(row)
	return None, 'figure not written'


def _stage(sample, pattern):
	# Substitute only {sample}; figure patterns keep their {ext} placeholder so
	# each consumer can pick png (HTML) or pdf (slides) from the same entry.
	return pattern.replace('{sample}', sample)


def cellranger_metrics(samples, sampledir):
	"""Cell Ranger metrics_summary.csv per sample, concatenated.

	Missing files are reported, not silently skipped -- an absent metrics table
	usually means the Cell Ranger path is wrong.
	"""
	frames, missing = [], []
	for sid, rel in samples['cellranger'].items():
		path = os.path.join(sampledir, rel, 'metrics_summary.csv')
		if not os.path.exists(path):
			missing.append(sid)
			continue
		df = pd.read_csv(path)
		df.insert(0, 'sampleid', sid)
		frames.append(df)
	table = pd.concat(frames, ignore_index=True) if frames else None
	return table, missing


def collect(samples, sampledir, config, nf_samples, callers, status_file=None):
	"""Every table and figure path both reports need.

	`status_file` is result/qc_status.csv. Without it (a caller from before
	v0.3.4) the reports are built as if every step had succeeded.
	"""
	ids = samples['sample'].tolist()
	metrics, metrics_missing = cellranger_metrics(samples, sampledir)
	status = read_status(status_file)

	data = {
		'samples': samples,
		'sample_ids': ids,
		'config': config,
		'seed': config.get('seed'),
		'nf_samples': list(nf_samples),
		'nf_missing': [s for s in ids if s not in set(nf_samples)],
		'callers': list(callers),
		'decider': config['doublet']['decider'],
		'ambient_method': config['ambient']['method'],
		'ambient_compare': list(config['ambient'].get('compare', [])),
		'cellranger_metrics': metrics,
		'cellranger_metrics_missing': metrics_missing,
		'status': status,
		'ambient': _every_sample(_read_concat([_stage(s, 'ambient/{sample}_contamination.txt') for s in ids]), ids),
		'barcoderank': _every_sample(_read_concat([_stage(s, 'barcoderank/{sample}_knee.txt') for s in ids]), ids),
		'filter_ncell': _every_sample(_read_concat([_stage(s, 'filterbycount/{sample}_filter_ncell.txt') for s in ids]), ids),
		'doublet_summary': _every_sample(_read_concat([_stage(s, 'filterdoublet/{sample}_doublet_summary.txt') for s in ids]), ids),
		'concordance': _every_sample(_read_concat([_stage(s, 'filterdoublet/{sample}_doublet_concordance.txt') for s in ids]), ids),
		# Each caller's own count, read from the caller rather than from
		# filterdoublet: a caller that ran on a sample whose decider then failed
		# still has a number, and the metrics table should carry it.
		'doublet_ratio': _read_concat([_stage(s, f'{c}/{{sample}}_doublet_ratio.txt') for c in callers for s in ids]),
	}
	data.update(sample_outcome(status, ids))

	data['figures'] = {
		'barcoderank': {s: _stage(s, 'barcoderank/{sample}_barcoderank.{ext}') for s in ids},
		'ambient': {s: _stage(s, 'ambient/{sample}_ambient.{ext}') for s in ids},
		'violin_before': {s: _stage(s, 'filterbycount/{sample}_violin_before.{ext}') for s in ids},
		'violin_after': {s: _stage(s, 'filterbycount/{sample}_violin_after.{ext}') for s in ids},
		'nf': {s: _stage(s, 'nuclear_fraction/{sample}_nf_umi.{ext}') for s in nf_samples},
	}
	if 'doubletfinder' in callers:
		data['figures']['doubletfinder_pANN'] = {s: _stage(s, 'doubletfinder/{sample}_pANN.{ext}') for s in ids}
		data['figures']['doubletfinder_umap'] = {s: _stage(s, 'doubletfinder/{sample}_umap.{ext}') for s in ids}
	if 'scdblfinder' in callers:
		data['figures']['scdblfinder_score'] = {s: _stage(s, 'scdblfinder/{sample}_score.{ext}') for s in ids}

	data['nf_summary'] = _every_sample(nf_summary(nf_samples), [s for s in ids if s in set(nf_samples)])
	data['cascade'] = cascade(data, ids)
	data['metrics'] = metrics_table(data, ids)
	data['caveats'] = caveats(data)
	return data


def sample_outcome(status, ids):
	"""Who reached result/, who did not and why, and which steps fell back."""
	out = {'included': list(ids), 'excluded': {}, 'fallbacks': [], 'problems': None}
	if status is None:
		return out
	res = status[status['step'] == 'result'].set_index('sample')
	out['included'] = [s for s in ids if s in res.index and res.loc[s, 'status'] == qcstatus.OK]
	out['excluded'] = {
		s: res.loc[s, 'message'].removeprefix('excluded: ') for s in ids
		if s in res.index and res.loc[s, 'status'] != qcstatus.OK}
	fb = status[status['status'] == qcstatus.FALLBACK]
	out['fallbacks'] = [(r['sample'], r['step'], r['message']) for _, r in fb.iterrows()]
	# Everything that is not a plain `ok`, for the status section of both reports.
	bad = status[(status['status'] != qcstatus.OK) & (status['step'] != 'result')]
	out['problems'] = bad.reset_index(drop=True) if len(bad) else None
	return out


def nf_summary(nf_samples):
	"""Per-sample nuclear-fraction summary, for the metrics table.

	The step writes one row per barcode; the cohort table wants one number per
	sample. Samples without a Cell Ranger BAM have no file and are simply absent,
	which is why the metrics row shows blanks rather than zeros for them.
	"""
	rows = []
	for sid in nf_samples:
		path = _stage(sid, 'nuclear_fraction/{sample}.txt.gz')
		if not computed(path):
			continue
		df = pd.read_csv(path, sep='\t', header=0)
		nf = pd.to_numeric(df['nuclear_fraction'], errors='coerce')
		rows.append({
			'sampleid': sid,
			'ncell': int(len(df)),
			'median': float(np.nanmedian(nf)) if len(nf) else np.nan,
			'q25': float(np.nanquantile(nf, 0.25)) if len(nf) else np.nan,
			'q75': float(np.nanquantile(nf, 0.75)) if len(nf) else np.nan,
			'n_missing': int(nf.isna().sum()),
			})
	return pd.DataFrame(rows) if rows else None


def metrics_table(data, ids):
	"""Every scalar the run produced, one row per sample.

	Written to result/metrics.csv. The reports are for reading; this is for
	joining -- a cohort summary, a spreadsheet, or a downstream script that
	should not have to parse six stage-specific stats files, and must never have
	to scrape a number out of the HTML.

	Column names are namespaced by stage, and any name that depends on a
	configured method carries that method in it, so adding a caller or a backend
	adds columns instead of changing the meaning of existing ones.
	"""
	def sub(table, sid):
		if table is None or 'sampleid' not in getattr(table, 'columns', []):
			return None
		m = table[table['sampleid'] == sid]
		return m if len(m) else None

	rows = []
	for sid in ids:
		row = {'sampleid': sid}

		cr = sub(data['cellranger_metrics'], sid)
		if cr is not None:
			for col in cr.columns:
				if col != 'sampleid':
					row['cellranger_' + col.strip().lower().replace(' ', '_')] = cr[col].iloc[0]

		knee = sub(data['barcoderank'], sid)
		if knee is not None:
			for col in knee.columns:
				if col != 'sampleid':
					row['barcoderank_' + col] = knee[col].iloc[0]

		amb = sub(data['ambient'], sid)
		if amb is not None:
			for _, r in amb.iterrows():
				name = r.get('method', 'ambient')
				for col in amb.columns:
					if col not in ('sampleid', 'method'):
						row[f'ambient_{name}_{col}'] = r[col]

		fil = sub(data['filter_ncell'], sid)
		if fil is not None:
			for col in fil.columns:
				if col != 'sampleid':
					row['filter_' + col] = fil[col].iloc[0]

		dbl = sub(data['doublet_summary'], sid)
		if dbl is not None and 'caller' in dbl:
			for _, r in dbl.dropna(subset=['caller']).iterrows():
				caller = r['caller']
				row[f'doublet_{caller}_ndoublet'] = r.get('ndoublet')
				row[f'doublet_{caller}_frac'] = r.get('frac_doublet')
				if not pd.isna(r.get('is_decider')) and bool(r.get('is_decider')):
					row['doublet_decider'] = caller
					row['doublet_ncell_before'] = r.get('ncell_before')
					row['doublet_ncell_after'] = r.get('ncell_after')
		ratio = sub(data.get('doublet_ratio'), sid)
		if ratio is not None:
			for _, r in ratio.iterrows():
				caller = r['caller']
				if pd.isna(row.get(f'doublet_{caller}_ndoublet', np.nan)):
					row[f'doublet_{caller}_ndoublet'] = r.get('ndoublet')
					n = r.get('ncell_before')
					row[f'doublet_{caller}_frac'] = r.get('ndoublet') / n if n else np.nan

		conc = sub(data['concordance'], sid)
		if conc is not None and 'caller_a' in conc:
			for _, r in conc.dropna(subset=['caller_a']).iterrows():
				pair = f"{r['caller_a']}_vs_{r['caller_b']}"
				row[f'concordance_{pair}_kappa'] = r.get('kappa')
				row[f'concordance_{pair}_both_doublet'] = r.get('both_doublet')

		nf = sub(data.get('nf_summary'), sid)
		if nf is not None:
			for col in nf.columns:
				if col != 'sampleid':
					row['nf_' + col] = nf[col].iloc[0]

		casc = sub(data['cascade'], sid)
		if casc is not None and 'frac_retained' in casc:
			row['frac_retained'] = casc['frac_retained'].iloc[0]

		# Last, so the columns of a run with no failures keep their v0.3.3
		# positions. One status column per step, plus whether the sample reached
		# result/ and, if not, why.
		st = data.get('status')
		if st is not None:
			for _, r in st[st['sample'] == sid].iterrows():
				if r['step'] != 'result':
					row['status_' + r['step'].replace('.', '_')] = r['status']
		row['qc_included'] = sid in data['included']
		row['qc_excluded_reason'] = data['excluded'].get(sid, '')

		rows.append(row)
	table = pd.DataFrame(rows)
	# A count column with an NA in it would otherwise become float and print a
	# cell count as 512.0. Columns whose every value is an integer stay integer.
	def is_int(v):
		return isinstance(v, (int, np.integer)) and not isinstance(v, (bool, np.bool_))
	for c in table.columns:
		vals = [r[c] for r in rows if c in r and not pd.isna(r[c])]
		if vals and all(is_int(v) for v in vals):
			table[c] = table[c].astype('Int64')
	return table


def cascade(data, ids):
	"""Cells surviving each stage, one row per sample.

	This is the table that makes the pipeline auditable: every exclusion appears
	as a difference between two adjacent columns.
	"""
	def count(m, col):
		if not len(m) or col not in m:
			return pd.NA
		v = m[col].iloc[0]
		return pd.NA if pd.isna(v) else int(v)

	rows = []
	bcr = data['barcoderank']
	fil = data['filter_ncell']
	dbl = data['doublet_summary']
	for sid in ids:
		row = {'sampleid': sid}
		if bcr is not None:
			row['cellranger_cells'] = count(bcr[bcr['sampleid'] == sid], 'n_called_cells')
		if fil is not None:
			m = fil[fil['sampleid'] == sid]
			row['after_filterbycount'] = count(m, 'ncell_after')
			row['removed_by_count'] = count(m, 'ncell_removed')
		if dbl is not None:
			m = dbl[(dbl['sampleid'] == sid)]
			dec = m[m['is_decider'].fillna(False).astype(bool)] if 'is_decider' in m else m
			dec = dec if len(dec) else m.iloc[0:0]
			row['after_doublet'] = count(dec, 'ncell_after')
			row['removed_by_doublet'] = count(dec, 'ndoublet')
		cells, after = row.get('cellranger_cells', pd.NA), row.get('after_doublet', pd.NA)
		row['frac_retained'] = round(after / cells, 4) if not (pd.isna(cells) or pd.isna(after)) and cells else np.nan
		row['in_result'] = 'yes' if sid in data['included'] else 'no'
		rows.append(row)
	table = pd.DataFrame(rows)
	for c in ('cellranger_cells', 'after_filterbycount', 'removed_by_count', 'after_doublet', 'removed_by_doublet'):
		if c in table:
			table[c] = table[c].astype('Int64')
	return table


def caveats(data):
	"""Limitations that must travel with the results.

	Stated in the output rather than left in the source, because a report that
	presents QC numbers without them invites over-reading.
	"""
	items = []
	if data['excluded']:
		items.append(
			f"{len(data['excluded'])} of {len(data['sample_ids'])} sample(s) failed a required "
			'step and are NOT in result/: '
			+ '; '.join(f'{s} ({why})' for s, why in data['excluded'].items())
			+ '. They still appear in every table, with NA for what was not computed; '
			'result/qc_status.csv has every step of every sample.'
			)
	if data['fallbacks']:
		items.append(
			'Fallbacks were used, so these samples were not processed exactly as '
			'configured: ' + '; '.join(f'{s} {step}: {msg}' for s, step, msg in data['fallbacks']) + '.'
			)
	items += [
		'Homotypic doublets are not modelled (modelHomotypic is deliberately not '
		'called), so the expected-doublet count over-estimates the DETECTABLE '
		'doublet count and the doublet step removes slightly more cells than the '
		'true heterotypic count. The bias direction is known and constant.',
		'The expected doublet rate is the 10x rule of thumb '
		f"(rate={data['config']['doublet'].get('rate')}, "
		f"capacity={data['config']['doublet'].get('capacity')} cells per reaction), "
		'not a measurement for this library.',
		'Cell calling is Cell Ranger EmptyDrops. CellQC does not re-call cells; the '
		'barcode rank plot is diagnostic only.',
	]
	if data['nf_samples']:
		items.append(
			'The nuclear fraction is reported but NOT used for filtering. '
			'DropletQC-style empty-drop and damaged-cell thresholds are sample- and '
			'tissue-dependent, so applying them automatically across a cohort would '
			'be unreviewed auto-filtering.'
			)
	if data['nf_missing']:
		items.append(
			f"No Cell Ranger BAM for {len(data['nf_missing'])} sample(s) "
			f"({', '.join(data['nf_missing'])}), so the nuclear fraction was not "
			'computed for them. A missing column is not a zero.'
			)
	if len(data['callers']) > 1:
		items.append(
			f"Two doublet callers ran; only {data['decider']} removed cells. The "
			'concordance table is a consistency measure, NOT evidence that either '
			'caller is correct -- with no ground-truth doublets neither can be shown '
			'superior on this data.'
			)
	if data['ambient_compare']:
		items.append(
			f"Ambient correction applied: {data['ambient_method']}. "
			f"Also estimated for comparison: {', '.join(data['ambient_compare'])} "
			'(these did NOT modify counts).'
			)
	if data['cellranger_metrics_missing']:
		items.append(
			'metrics_summary.csv missing for: '
			f"{', '.join(data['cellranger_metrics_missing'])}."
			)
	items.append(
		f"Random seed {data['seed']} was set for every stochastic step. "
		'v0.1.0 set no seed and its doublet calls were not reproducible.'
		)
	return items
