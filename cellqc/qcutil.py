# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:
"""Shared helpers for cellqc Snakemake scripts.

Lives inside the installed package (not next to the scripts) because Snakemake
copies `script:` targets into .snakemake/scripts/ before running them, which
breaks sibling imports.
"""

import os
import random

import numpy as np

# --- Figure policy -----------------------------------------------------------
# Every figure is emitted twice from one canvas: a vector PDF whose text stays
# editable (Type 42 / TrueType, never Type 3 outlines) for slides and figures,
# and a raster PNG for the HTML report, which cannot inline a PDF in an <img>.
# Dense layers (scatter of ~1e4 points) are rasterized inside the PDF at 500 dpi
# so the file stays small without turning the axis labels into pixels.

RASTER_DPI = 500
PNG_DPI = 200


def setup_matplotlib():
	"""Configure matplotlib for editable-text vector output. Call before pyplot."""
	import matplotlib

	matplotlib.use('Agg')
	matplotlib.rcParams.update({
		'pdf.fonttype': 42,      # TrueType: text stays selectable/editable
		'ps.fonttype': 42,
		'svg.fonttype': 'none',
		'figure.dpi': 100,
		'savefig.bbox': 'tight',
		'axes.spines.top': False,
		'axes.spines.right': False,
		'font.size': 9,
	})
	return matplotlib


def savefig(fig, stem, png=True):
	"""Save `fig` as <stem>.pdf (vector, editable text) and optionally <stem>.png.

	`stem` is a path without extension. Returns the list of paths written.
	"""
	written = []
	os.makedirs(os.path.dirname(os.path.abspath(stem)) or '.', exist_ok=True)
	fig.savefig(f'{stem}.pdf', dpi=RASTER_DPI)
	written.append(f'{stem}.pdf')
	if png:
		fig.savefig(f'{stem}.png', dpi=PNG_DPI)
		written.append(f'{stem}.png')
	return written


# --- Reproducibility ---------------------------------------------------------

def set_seed(seed):
	"""Seed every RNG a Python QC step can reach.

	v0.1.0 seeded nothing, so its stochastic steps were not reproducible. Any
	script that clusters, embeds, or subsamples must call this.
	"""
	seed = int(seed)
	random.seed(seed)
	np.random.seed(seed)
	os.environ['PYTHONHASHSEED'] = str(seed)
	return seed


# --- Cell Ranger layout ------------------------------------------------------

BAM_NAME = 'possorted_genome_bam.bam'


def has_bam(cellranger_dir):
	"""True when a Cell Ranger dir carries an indexed BAM.

	Drives the automatic skip of the nuclear-fraction step: no BAM means the
	statistic cannot be computed, so the step is dropped from the DAG for that
	sample rather than failing the run.
	"""
	bam = os.path.join(cellranger_dir, BAM_NAME)
	return os.path.exists(bam) and os.path.exists(bam + '.bai')


def raw_h5(cellranger_dir):
	return os.path.join(cellranger_dir, 'raw_feature_bc_matrix.h5')


def filtered_h5(cellranger_dir):
	return os.path.join(cellranger_dir, 'filtered_feature_bc_matrix.h5')


def metrics_csv(cellranger_dir):
	return os.path.join(cellranger_dir, 'metrics_summary.csv')


# --- .obs column names -------------------------------------------------------
# Namespaced per caller so downstream code never has to guess which method
# produced a column, and so swapping the deciding caller changes no schema.

OBS_DF_SCORE = 'doubletfinder_pANN'
OBS_DF_CLASS = 'doubletfinder_class'
OBS_SDB_SCORE = 'scdblfinder_score'
OBS_SDB_CLASS = 'scdblfinder_class'
OBS_NF = 'nuclear_fraction'

DOUBLET_CALLERS = ('doubletfinder', 'scdblfinder')

CALLER_COLUMNS = {
	'doubletfinder': (OBS_DF_SCORE, OBS_DF_CLASS),
	'scdblfinder': (OBS_SDB_SCORE, OBS_SDB_CLASS),
}


def expected_doublet_rate(ncell, nreaction, rate, capacity):
	"""10x multiplet-rate rule of thumb, as used since v0.1.0.

	rate * ncell / (nreaction * capacity) -- with the v0.1.0 defaults
	(rate=0.1, capacity=13000) this is the familiar ~0.8% per 1,000 cells
	recovered. The two constants were hard-coded in v0.1.0; they are config
	parameters now so the assumption is visible, but the defaults reproduce
	v0.1.0 exactly.

	Returns (ratio, n_expected). Not adjusted for homotypic doublets -- see the
	caveat carried in the reports.
	"""
	ratio = round(rate * ncell / (nreaction * capacity), 2)
	return ratio, int(round(ratio * ncell))
