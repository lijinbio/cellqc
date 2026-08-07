# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:
"""Assemble everything both reports show, from the files the workflow wrote.

The HTML report and the PDF slide deck call this and nothing else, so the two
cannot drift into describing different runs. Nothing here recomputes an
analysis: every number is read back from the stage that produced it.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd


def _read_concat(paths):
	frames = [pd.read_csv(p, sep='\t', header=0) for p in paths if os.path.exists(p)]
	frames = [f for f in frames if len(f)]
	return pd.concat(frames, ignore_index=True) if frames else None


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


def collect(samples, sampledir, config, nf_samples, callers):
	"""Every table and figure path both reports need."""
	ids = samples['sample'].tolist()
	metrics, metrics_missing = cellranger_metrics(samples, sampledir)

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
		'ambient': _read_concat([_stage(s, 'ambient/{sample}_contamination.txt') for s in ids]),
		'barcoderank': _read_concat([_stage(s, 'barcoderank/{sample}_knee.txt') for s in ids]),
		'filter_ncell': _read_concat([_stage(s, 'filterbycount/{sample}_filter_ncell.txt') for s in ids]),
		'doublet_summary': _read_concat([_stage(s, 'result/{sample}_doublet_summary.txt') for s in ids]),
		'concordance': _read_concat([_stage(s, 'result/{sample}_doublet_concordance.txt') for s in ids]),
	}

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

	data['nf_summary'] = nf_summary(nf_samples)
	data['cascade'] = cascade(data, ids)
	data['metrics'] = metrics_table(data, ids)
	data['caveats'] = caveats(data)
	return data


def nf_summary(nf_samples):
	"""Per-sample nuclear-fraction summary, for the metrics table.

	The step writes one row per barcode; the cohort table wants one number per
	sample. Samples without a Cell Ranger BAM have no file and are simply absent,
	which is why the metrics row shows blanks rather than zeros for them.
	"""
	rows = []
	for sid in nf_samples:
		path = _stage(sid, 'nuclear_fraction/{sample}.txt.gz')
		if not os.path.exists(path):
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
		if dbl is not None:
			for _, r in dbl.iterrows():
				caller = r['caller']
				row[f'doublet_{caller}_ndoublet'] = r.get('ndoublet')
				row[f'doublet_{caller}_frac'] = r.get('frac_doublet')
				if r.get('is_decider'):
					row['doublet_decider'] = caller
					row['doublet_ncell_before'] = r.get('ncell_before')
					row['doublet_ncell_after'] = r.get('ncell_after')

		conc = sub(data['concordance'], sid)
		if conc is not None:
			for _, r in conc.iterrows():
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

		rows.append(row)
	return pd.DataFrame(rows)


def cascade(data, ids):
	"""Cells surviving each stage, one row per sample.

	This is the table that makes the pipeline auditable: every exclusion appears
	as a difference between two adjacent columns.
	"""
	rows = []
	bcr = data['barcoderank']
	fil = data['filter_ncell']
	dbl = data['doublet_summary']
	for sid in ids:
		row = {'sampleid': sid}
		if bcr is not None:
			m = bcr[bcr['sampleid'] == sid]
			row['cellranger_cells'] = int(m['n_called_cells'].iloc[0]) if len(m) else np.nan
		if fil is not None:
			m = fil[fil['sampleid'] == sid]
			if len(m):
				row['after_filterbycount'] = int(m['ncell_after'].iloc[0])
				row['removed_by_count'] = int(m['ncell_removed'].iloc[0])
		if dbl is not None:
			m = dbl[(dbl['sampleid'] == sid)]
			if len(m):
				dec = m[m.get('is_decider', False) == True] if 'is_decider' in m else m
				dec = dec if len(dec) else m
				row['after_doublet'] = int(dec['ncell_after'].iloc[0]) if 'ncell_after' in dec else np.nan
				row['removed_by_doublet'] = int(dec['ndoublet'].iloc[0]) if 'ndoublet' in dec else np.nan
		if row.get('cellranger_cells') and row.get('after_doublet'):
			row['frac_retained'] = round(row['after_doublet'] / row['cellranger_cells'], 4)
		rows.append(row)
	return pd.DataFrame(rows)


def caveats(data):
	"""Limitations that must travel with the results.

	Stated in the output rather than left in the source, because a report that
	presents QC numbers without them invites over-reading.
	"""
	items = [
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
