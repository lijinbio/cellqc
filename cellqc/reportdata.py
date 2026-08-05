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
		'decider': config['doublet'].get('decider') if not config['doublet'].get('skip') else None,
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

	data['cascade'] = cascade(data, ids)
	data['caveats'] = caveats(data)
	return data


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
		'Cell calling is Cell Ranger EmptyDrops. cellqc does not re-call cells; the '
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
