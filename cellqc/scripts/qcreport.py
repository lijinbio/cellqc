# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:
"""Self-contained HTML QC report.

Every figure is inlined as a base64 data URI so the file can be emailed or
archived on its own. PNG twins are used here because an <img> cannot render the
vector PDFs; the PDFs go into the slide deck.
"""

import base64
import mimetypes
import os
import sys
from pathlib import Path

import pandas as pd
from jinja2 import Environment, FileSystemLoader

from cellqc import __version__, reportdata

samples = snakemake.params['samples']
sampledir = snakemake.params['sampledir']
config = snakemake.params['config']
nf_samples = snakemake.params['nf_samples']
callers = snakemake.params['callers']
outfile = snakemake.output['html']
metricsfile = snakemake.output['metrics']

TEMPLATE_DIR = Path(__file__).parent / 'template'


# From snakemake/report/__init__.py
def data_uri(data, filename, encoding='utf8', mime='text/plain'):
	data = base64.b64encode(data)
	return 'data:{mime};charset={charset};filename={filename};base64,{data}'.format(
		filename=filename, mime=mime, charset=encoding, data=data.decode('utf-8'))


def data_uri_from_file(file, defaultenc='utf8'):
	file = str(file)
	mime, encoding = mimetypes.guess_type(file)
	if mime is None:
		mime = 'text/plain'
		print(f'Could not detect mimetype for {file}, assuming text/plain.', file=sys.stderr)
	if encoding is None:
		encoding = defaultenc
	with open(file, 'rb') as f:
		return data_uri(f.read(), os.path.basename(file), encoding, mime)


def get_resource_as_string(path):
	return open(TEMPLATE_DIR / path).read()


def to_html(df):
	if df is None or not len(df):
		return None
	return df.to_html(index=False, na_rep='', float_format=lambda v: f'{v:,.4g}')


def figure_sections(data):
	"""Ordered figure blocks, each a list of (sample, data_uri)."""
	spec = [
		('barcoderank', 'Barcode rank',
			'Cell Ranger EmptyDrops call (dashed) against the UMI-rank curve. '
			'Diagnostic only — CellQC does not re-call cells.'),
		('ambient', f"Ambient RNA — {data['ambient_method']} (applied)",
			'Estimated contamination fraction. This correction modified the counts.'),
		('violin_before', 'QC metrics before filtering',
			'Dashed lines mark the applied thresholds.'),
		('violin_after', 'QC metrics after filtering', ''),
		('nf', 'Nuclear fraction vs UMI depth',
			'Intronic / (intronic + exonic) reads per cell, coloured by % mitochondrial '
			'UMI before ambient correction where the reference has mitochondrial genes. '
			'Reported only — NOT used for filtering.'),
		('doubletfinder_pANN', 'DoubletFinder — pANN', ''),
		('doubletfinder_umap', 'DoubletFinder — calls on UMAP', ''),
		('scdblfinder_score', 'scDblFinder — score',
			'Run with the same expected doublet rate as DoubletFinder. Diagnostic only.'),
		]
	out = []
	for key, title, blurb in spec:
		paths = data['figures'].get(key, {})
		imgs = []
		for sid, pattern in paths.items():
			png = pattern.format(ext='png')
			if os.path.exists(png):
				imgs.append((sid, data_uri_from_file(png)))
		if imgs:
			out.append({'title': title, 'blurb': blurb, 'images': imgs})
	return out


def main():
	data = reportdata.collect(samples, sampledir, config, nf_samples, callers)

	env = Environment(loader=FileSystemLoader(str(TEMPLATE_DIR)))
	env.filters['get_resource_as_string'] = get_resource_as_string
	template = env.get_template('index.html.jinja2')

	html = template.render(
		version=__version__,
		seed=data['seed'],
		nsample=len(data['sample_ids']),
		ambient_method=data['ambient_method'],
		ambient_compare=data['ambient_compare'],
		callers=data['callers'],
		decider=data['decider'],
		cascade=to_html(data['cascade']),
		cellrangersummary=to_html(data['cellranger_metrics']),
		barcoderank=to_html(data['barcoderank']),
		ambient=to_html(data['ambient']),
		filterncell=to_html(data['filter_ncell']),
		doubletsummary=to_html(data['doublet_summary']),
		concordance=to_html(data['concordance']),
		sections=figure_sections(data),
		caveats=data['caveats'],
		)

	with open(outfile, mode='w', encoding='utf-8') as f:
		f.write(html)
	print(f'[qcreport] wrote {outfile}', flush=True)

	# The machine-readable twin of the report: every scalar the run produced,
	# one row per sample. Written here so it is assembled from the same
	# reportdata.collect() call the HTML is, and cannot disagree with it.
	data['metrics'].to_csv(metricsfile, index=False)
	print(f'[qcreport] wrote {metricsfile} '
		f"({data['metrics'].shape[0]} samples x {data['metrics'].shape[1]} metrics)", flush=True)


if __name__ == '__main__':
	main()
