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

# Cell Ranger's metrics table has ~20 columns and will not fit a 16:9 slide as
# one row. It is transposed and folded into two metric/value pairs side by side,
# which fills the slide at a readable size instead of leaving a strip of 5pt
# type in the middle of an empty frame.

# Headers straight off the stats files are snake_case and long enough to push the
# cascade table off the right-hand edge of the slide. The table is the one a
# reader is meant to audit, so it gets real column names.
CASCADE_HEADERS = {
	'sampleid': 'Sample',
	'cellranger_cells': 'Cell Ranger cells',
	'after_filterbycount': 'After count filter',
	'removed_by_count': 'Removed (count)',
	'after_doublet': 'After doublet',
	'removed_by_doublet': 'Removed (doublet)',
	'frac_retained': 'Retained',
	}


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


def tex_table(df, max_rows=None, headers=None, align=None):
	"""A booktabs tabular. Returns None for an empty frame so the template can
	drop the frame rather than emit an empty box.

	`headers` renames columns for display; `align` is a per-column alignment
	string (defaults to left for the first column and right for the rest, which
	is what makes columns of counts readable).
	"""
	if df is None or not len(df):
		return None
	df = df.copy()
	if max_rows:
		df = df.head(max_rows)
	cols = list(df.columns)
	labels = [(headers or {}).get(c, c) for c in cols]
	if align is None:
		align = 'l' + 'r' * (len(cols) - 1)
	lines = [r'\begin{tabular}{' + align + '}', r'\toprule']
	lines.append(' & '.join(rf'\textbf{{{tex_escape(c)}}}' for c in labels) + r' \\')
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


def metrics_table(row):
	"""Cell Ranger metrics for one sample as a two-up metric/value tabular.

	Transposed (20 columns do not fit a slide) and folded in half so the frame
	is filled at a legible size. `sampleid` is dropped: it is the frame title.
	"""
	cols = [c for c in row.columns if c != 'sampleid']
	pairs = [(c, row.iloc[0][c]) for c in cols]
	half = (len(pairs) + 1) // 2
	left, right = pairs[:half], pairs[half:]
	right += [('', '')] * (len(left) - len(right))
	folded = pd.DataFrame({
		'Metric': [k for k, _ in left],
		'Value': [v for _, v in left],
		'Metric ': [k for k, _ in right],
		'Value ': [v for _, v in right],
		})
	return tex_table(folded, align='lrlr')


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
		m = None
		if metrics is not None:
			row = metrics[metrics['sampleid'] == sid]
			if len(row):
				m = metrics_table(row)
		per_sample.append({
			'label': tex_escape(sid),
			'metrics': m,
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
		cascade_table=tex_table(data['cascade'], headers=CASCADE_HEADERS),
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
