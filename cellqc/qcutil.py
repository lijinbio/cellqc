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
# Dense layers (scatter of ~1e4 points) are rasterized inside the PDF at 600 dpi
# so the file stays small without turning the axis labels into pixels.
#
# PNG_DPI is 300 -- print resolution. The PNGs are base64-inlined into the HTML
# report, so this is the one figure setting that trades report size for
# legibility; 200 dpi was visibly soft when a reader zoomed into a violin plot.

RASTER_DPI = 600
PNG_DPI = 300


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


# --- Mitochondrial content ---------------------------------------------------
# One definition, used by `filterbycount` (via scanpy's qc_vars=['mt']) and by
# the nuclear-fraction scatter, so the number a cell is filtered on and the
# number colouring that cell in the plot cannot drift apart.

MITO_PREFIXES = ('MT-', 'mt-')


def mito_percent(adata):
	"""Percent of a cell's UMI in mitochondrial genes.

	Same gene pattern and same formula as scanpy's `pct_counts_mt`, computed
	directly so a caller that only needs this one metric does not have to run
	`calculate_qc_metrics` over the whole matrix. Returns `(pct, n_mito_genes)`;
	`pct` is NaN for a barcode with no counts, and all-zero when nothing matched
	the pattern -- `n_mito_genes` is what tells those two cases apart.
	"""
	mt = np.fromiter(
		(str(n).startswith(MITO_PREFIXES) for n in adata.var_names),
		dtype=bool, count=adata.n_vars,
		)
	total = np.asarray(adata.X.sum(axis=1)).ravel().astype(float)
	mt_counts = np.asarray(adata.X[:, mt].sum(axis=1)).ravel().astype(float)
	with np.errstate(invalid='ignore', divide='ignore'):
		pct = np.where(total > 0, 100.0 * mt_counts / total, np.nan)
	return pct, int(mt.sum())


# --- Matrix annotations ------------------------------------------------------

OBS_INDEX_NAME = 'barcode'
VAR_INDEX_NAME = 'gene'


def write_obs_var(adata, obs_file, var_file):
	"""Write `.obs` and `.var` beside the matrix as gzipped TSVs.

	Every matrix cellqc keeps gets these, so the per-cell QC metrics and the
	feature table can be read (and joined on `barcode` / `gene`) without opening
	the `.h5ad` at all -- by R, by a spreadsheet, or by someone who does not have
	anndata installed. The index is named rather than left blank, because an
	unnamed first column is what makes a table like this ambiguous to reload.
	"""
	obs = adata.obs.copy()
	obs.index = obs.index.astype(str)
	obs.index.name = OBS_INDEX_NAME
	obs.to_csv(obs_file, sep='\t', index=True, compression='gzip')

	var = adata.var.copy()
	var.index = var.index.astype(str)
	var.index.name = VAR_INDEX_NAME
	var.to_csv(var_file, sep='\t', index=True, compression='gzip')
	return obs_file, var_file


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

	rate * ncell / (nreaction * capacity) -- linear in the number of cells
	recovered. The linearity is not an approximation of convenience: cells are
	loaded at limiting dilution, so occupancy is Poisson with mean lambda, and
	the fraction of occupied droplets holding two or more cells is
	1 - lambda/(exp(lambda) - 1), which is lambda/2 to first order. lambda is
	proportional to the cells loaded, hence the multiplet fraction is
	proportional to the cells recovered, and the rate is quoted per thousand
	cells rather than as one number.

	With the v0.1.0 defaults (rate=0.1, capacity=13000) this is 0.77% per 1,000
	cells recovered, i.e. the ~0.8% per 1,000 that 10x publishes. scDblFinder's
	default assumes ~1% per 1,000 (rate=0.1, capacity=10000 here). The constants
	were hard-coded in v0.1.0; they are config parameters now so the assumption
	is visible, and the defaults reproduce v0.1.0 exactly.

	`nreaction` divides the fraction because pooled reactions are separate
	emulsions: a cell from one reaction cannot share a droplet with a cell from
	another.

	Two known biases, both upward: the linear form is the small-lambda limit and
	sits slightly above the exact Poisson expression at high yields, and
	homotypic doublets are not modelled -- see the caveat carried in the reports.

	References: Bloom (2018) PeerJ 6:e5578 (exact Poisson treatment);
	10x Genomics Chromium user guides (multiplet rate vs cell recovery);
	Germain et al. (2021) F1000Research 10:979 (scDblFinder's dbr).

	Returns (ratio, n_expected).
	"""
	ratio = round(rate * ncell / (nreaction * capacity), 2)
	return ratio, int(round(ratio * ncell))
