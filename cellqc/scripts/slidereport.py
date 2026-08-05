# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:
"""Presentation-ready PDF slide deck, rendered from a Jinja2 LaTeX template and
compiled with tectonic.

The template loops over the sample table and over a declared figure list, so a
new sample needs no template edit and a missing optional figure drops its frame
instead of failing the build.

Figures are embedded as the vector PDFs the analysis steps wrote, so text inside
the figures stays selectable and editable in the deck.
"""

import datetime
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from jinja2 import Environment, FileSystemLoader

from cellqc import __version__, reportdata

samples = snakemake.params['samples']
sampledir = snakemake.params['sampledir']
config = snakemake.params['config']
nf_samples = snakemake.params['nf_samples']
callers = snakemake.params['callers']
nowtimestr = snakemake.params['nowtimestr']
outfile = snakemake.output[0]

TEMPLATE_DIR = Path(__file__).parent / 'template'

# Cell Ranger's metrics table has 20 columns and will not fit a 16:9 slide.
# Split it rather than shrinking it to unreadable type.
METRICS_SPLIT = 10


def tex_escape(text):
	if text is None:
		return ''
	out = str(text)
	for a, b in (
		('\\', r'\textbackslash{}'), ('&', r'\&'), ('%', r'\%'), ('$', r'\$'),
		('#', r'\#'), ('_', r'\_'), ('{', r'\{'), ('}', r'\}'),
		('~', r'\textasciitilde{}'), ('^', r'\textasciicircum{}'),
		):
		out = out.replace(a, b)
	return out


def tex_table(df, max_rows=None, transpose=False):
	"""A booktabs tabular. Returns None for an empty frame so the template can
	drop the frame rather than emit an empty box."""
	if df is None or not len(df):
		return None
	df = df.copy()
	if max_rows:
		df = df.head(max_rows)
	if transpose:
		df = df.T.reset_index()
		df.columns = ['metric'] + [f'v{i}' for i in range(1, df.shape[1])]
	cols = list(df.columns)
	lines = [r'\begin{tabular}{' + 'l' * len(cols) + '}', r'\toprule']
	lines.append(' & '.join(rf'\textbf{{{tex_escape(c)}}}' for c in cols) + r' \\')
	lines.append(r'\midrule')
	for _, row in df.iterrows():
		vals = []
		for v in row:
			if isinstance(v, float):
				vals.append(tex_escape(f'{v:,.4g}'))
			else:
				vals.append(tex_escape(v))
		lines.append(' & '.join(vals) + r' \\')
	lines += [r'\bottomrule', r'\end{tabular}']
	return '\n'.join(lines)


def figure_entries(data, sid):
	"""Figures for one sample, in narrative order. Only those that exist."""
	figs = data['figures']
	spec = [
		('barcoderank', 'Barcode rank',
			'Cell Ranger EmptyDrops call (dashed line) against the UMI-rank curve. '
			'A clean library shows the call near the knee. Diagnostic only --- '
			'cellqc does not re-call cells.'),
		('ambient', f"Ambient RNA ({data['ambient_method']})",
			'Estimated contamination fraction. This correction WAS applied to the counts.'),
		('violin_before', 'QC metrics before filtering',
			'Dashed lines are the applied thresholds. Counts and genes on a log scale.'),
		('violin_after', 'QC metrics after filtering',
			'Same axes as the previous slide, after removing cells failing any threshold.'),
		('nf', 'Nuclear fraction vs UMI depth',
			'Intronic / (intronic + exonic) reads per cell against sequencing depth. '
			'High nuclear fraction at low depth suggests damaged cells or free nuclei. '
			'Reported only --- NOT used for filtering.'),
		('doubletfinder_pANN', 'DoubletFinder: pANN',
			'Proportion of artificial nearest neighbours, by call.'),
		('doubletfinder_umap', 'DoubletFinder: calls on UMAP',
			'Doublets should sit between transcriptionally distinct populations, '
			'not form their own island.'),
		('scdblfinder_score', 'scDblFinder: score',
			'Run with the same expected doublet rate as DoubletFinder so the two are '
			'compared under matched assumptions. Diagnostic only.'),
		]
	out = []
	for key, title, caption in spec:
		path = figs.get(key, {}).get(sid)
		if not path:
			continue
		pdf = path.format(ext='pdf')
		if os.path.exists(pdf):
			out.append({'title': title, 'caption': caption, 'path': os.path.abspath(pdf)})
	return out


def main():
	data = reportdata.collect(samples, sampledir, config, nf_samples, callers)

	env = Environment(
		loader=FileSystemLoader(str(TEMPLATE_DIR)),
		block_start_string='((*', block_end_string='*))',
		variable_start_string='(((', variable_end_string=')))',
		comment_start_string='((#', comment_end_string='#))',
		trim_blocks=True, lstrip_blocks=True, autoescape=False,
		)
	template = env.get_template('slides.tex.jinja2')

	metrics = data['cellranger_metrics']
	per_sample = []
	for sid in data['sample_ids']:
		m_a = m_b = None
		if metrics is not None:
			row = metrics[metrics['sampleid'] == sid]
			if len(row):
				cols = [c for c in row.columns if c != 'sampleid']
				first, second = cols[:METRICS_SPLIT], cols[METRICS_SPLIT:]
				m_a = tex_table(row[['sampleid'] + first], transpose=True)
				if second:
					m_b = tex_table(row[second], transpose=True)
		per_sample.append({
			'label': tex_escape(sid),
			'metrics_a': m_a,
			'metrics_b': m_b,
			'figures': figure_entries(data, sid),
			})

	rendered = template.render(
		title=tex_escape(f'cellqc QC report'),
		date=datetime.date.today().isoformat(),
		version=tex_escape(__version__),
		seed=data['seed'],
		nsample=len(data['sample_ids']),
		nf_ran=len(data['nf_samples']),
		ambient_method=tex_escape(data['ambient_method']),
		ambient_compare=[tex_escape(m) for m in data['ambient_compare']],
		callers=[tex_escape(c) for c in data['callers']],
		decider=tex_escape(data['decider']) if data['decider'] else None,
		cfg_mincount=config['filterbycount']['mincount'],
		cfg_minfeature=config['filterbycount']['minfeature'],
		cfg_mito=config['filterbycount']['mito'],
		cfg_rate=config['doublet']['rate'],
		cfg_capacity=config['doublet']['capacity'],
		cascade_table=tex_table(data['cascade']),
		samples=per_sample,
		caveats=[tex_escape(c) for c in data['caveats']],
		)

	outpath = Path(outfile).resolve()
	outpath.parent.mkdir(parents=True, exist_ok=True)
	texfile = outpath.with_suffix('.tex')
	texfile.write_text(rendered, encoding='utf-8')
	print(f'[slidereport] wrote {texfile}', flush=True)

	if shutil.which('tectonic') is None:
		raise RuntimeError(
			'tectonic not found on PATH. It is in envs/cellqc.yaml; activate the '
			'cellqc environment or install it to build the slide deck.'
			)

	with tempfile.TemporaryDirectory() as tmp:
		proc = subprocess.run(
			['tectonic', '--keep-logs', '--outdir', tmp, str(texfile)],
			capture_output=True, text=True,
			)
		if proc.returncode != 0:
			log = Path(tmp) / (texfile.stem + '.log')
			tail = log.read_text(errors='replace')[-3000:] if log.exists() else ''
			raise RuntimeError(
				f'tectonic failed (exit {proc.returncode}).\n'
				f'--- stderr ---\n{proc.stderr[-3000:]}\n--- log tail ---\n{tail}'
				)
		shutil.move(os.path.join(tmp, texfile.stem + '.pdf'), str(outpath))

	print(f'[slidereport] built {outpath}', flush=True)


if __name__ == '__main__':
	main()
