# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:
"""Nuclear fraction per cell barcode, from the Cell Ranger BAM.

	nuclear_fraction = n_intronic / (n_intronic + n_exonic)

This reimplements DropletQC::nuclear_fraction_tags(), whose only dependency was
a GitHub-only R package. Two deliberate differences from that reference, both
measured rather than assumed (see docs/v0.2.0_plan.md 4.2):

  1. DropletQC splits the genome into 100 tiles and queries each with
     scanBam(which=), so a read spanning a tile boundary is counted twice. This
     makes one pass per contig and counts every record once.
  2. DropletQC's default ScanBamParam applies no flag filter, so secondary and
     duplicate alignments are counted. That default is matched here so the
     comparison is like-for-like; counts under a primary-only filter are also
     recorded so the sensitivity of the statistic to that choice is visible
     instead of assumed away.

Intergenic reads (RE:A:I) are excluded from both numerator and denominator.
"""

import gzip
import multiprocessing as mp
import os

import numpy as np
import pandas as pd
import pysam
import scanpy as sc

from cellqc import qcutil

cellranger = snakemake.input['cellranger']
filtered_h5 = snakemake.input['filtered']
out_table, out_pdf, out_png = snakemake.output[0], snakemake.output[1], snakemake.output[2]
sampleid = snakemake.params['sampleid']
cbtag = snakemake.params['cbtag']
retag = snakemake.params['retag']
exontag = snakemake.params['exontag']
introntag = snakemake.params['introntag']
nthreads = max(1, int(snakemake.threads))

bam_path = os.path.join(cellranger, qcutil.BAM_NAME)


def cell_barcodes():
	"""The called-cell barcodes, in Cell Ranger's order."""
	adata = sc.read_10x_h5(filtered_h5)
	return list(adata.obs_names)


def count_contig(args):
	"""Count exonic/intronic reads per barcode on one contig.

	Returns (exon, intron, exon_primary, intron_primary) as dicts keyed by the
	barcode's index in `barcodes`.
	"""
	contig, bam_path, bc_index, cbtag, retag, exontag, introntag = args
	exon, intron = {}, {}
	exon_p, intron_p = {}, {}
	n_reads = n_tagged = 0
	with pysam.AlignmentFile(bam_path, 'rb') as bam:
		for read in bam.fetch(contig, multiple_iterators=True):
			n_reads += 1
			try:
				cb = read.get_tag(cbtag)
				re = read.get_tag(retag)
			except KeyError:
				continue
			idx = bc_index.get(cb)
			if idx is None:
				continue
			n_tagged += 1
			primary = not (read.is_secondary or read.is_supplementary)
			if re == exontag:
				exon[idx] = exon.get(idx, 0) + 1
				if primary:
					exon_p[idx] = exon_p.get(idx, 0) + 1
			elif re == introntag:
				intron[idx] = intron.get(idx, 0) + 1
				if primary:
					intron_p[idx] = intron_p.get(idx, 0) + 1
			# RE == intergenic is deliberately counted in neither
	return exon, intron, exon_p, intron_p, n_reads, n_tagged


def main():
	barcodes = cell_barcodes()
	bc_index = {bc: i for i, bc in enumerate(barcodes)}
	print(f'[nuclear_fraction] {sampleid}: {len(barcodes)} called-cell barcodes', flush=True)

	with pysam.AlignmentFile(bam_path, 'rb') as bam:
		contigs = [(c, bam.get_reference_length(c)) for c in bam.references]
	# Longest contig first so the slowest chunk starts earliest.
	contigs.sort(key=lambda x: -x[1])
	jobs = [
		(c, bam_path, bc_index, cbtag, retag, exontag, introntag)
		for c, _ in contigs
	]

	n = len(barcodes)
	exon = np.zeros(n, dtype=np.int64)
	intron = np.zeros(n, dtype=np.int64)
	exon_p = np.zeros(n, dtype=np.int64)
	intron_p = np.zeros(n, dtype=np.int64)
	total_reads = total_tagged = 0

	print(f'[nuclear_fraction] {sampleid}: scanning {len(jobs)} contigs on {nthreads} thread(s)', flush=True)
	if nthreads > 1:
		with mp.Pool(nthreads) as pool:
			results = pool.imap_unordered(count_contig, jobs)
			results = list(results)
	else:
		results = [count_contig(j) for j in jobs]

	for e, i, ep, ip, nr, nt in results:
		for d, arr in ((e, exon), (i, intron), (ep, exon_p), (ip, intron_p)):
			if d:
				idx = np.fromiter(d.keys(), dtype=np.int64, count=len(d))
				val = np.fromiter(d.values(), dtype=np.int64, count=len(d))
				np.add.at(arr, idx, val)
		total_reads += nr
		total_tagged += nt

	denom = exon + intron
	with np.errstate(invalid='ignore', divide='ignore'):
		nf = np.where(denom > 0, intron / denom, np.nan)
		denom_p = exon_p + intron_p
		nf_p = np.where(denom_p > 0, intron_p / denom_p, np.nan)

	# Never silently drop: report barcodes with no usable reads instead of
	# quietly emitting NaN.
	n_empty = int((denom == 0).sum())
	if n_empty:
		print(
			f'[nuclear_fraction] {sampleid}: WARNING {n_empty}/{n} barcodes have no '
			f'exonic or intronic reads; nuclear_fraction is NaN for these (kept, not dropped)',
			flush=True,
			)

	out = pd.DataFrame({
		'barcode': barcodes,
		'nuclear_fraction': nf,
		'exon': exon,
		'intron': intron,
		'nuclear_fraction_primary_only': nf_p,
		})
	with gzip.open(out_table, 'wt') as fh:
		out.to_csv(fh, sep='\t', index=False)

	frac_diff = np.nanmax(np.abs(nf - nf_p)) if n > 0 else np.nan
	print(
		f'[nuclear_fraction] {sampleid}: {total_reads} records scanned, '
		f'{total_tagged} assigned to called cells; '
		f'median NF={np.nanmedian(nf):.4f}; '
		f'max |NF - NF_primary_only|={frac_diff:.4f}',
		flush=True,
		)

	plot(barcodes, nf)


def plot(barcodes, nf):
	"""Nuclear fraction vs log10 total UMI, one point per called cell.

	Total UMI is taken from filtered_feature_bc_matrix.h5, i.e. the raw Cell
	Ranger counts before any ambient correction, so the axis means what it says.

	Descriptive only. cellqc does NOT filter on the nuclear fraction: the
	DropletQC empty-drop/damaged-cell thresholds are sample- and tissue-
	dependent, and applying them automatically across a cohort would be exactly
	the kind of unreviewed auto-filtering this pipeline avoids.
	"""
	matplotlib = qcutil.setup_matplotlib()
	import matplotlib.pyplot as plt

	adata = sc.read_10x_h5(filtered_h5)
	adata.var_names_make_unique()
	umi = np.asarray(adata.X.sum(axis=1)).ravel()
	order = pd.Index(adata.obs_names)
	umi = pd.Series(umi, index=order).reindex(barcodes).to_numpy()

	ok = np.isfinite(nf) & (umi > 0)
	fig, ax = plt.subplots(figsize=(4.6, 4.0))
	# Rasterize the dense scatter inside the vector PDF: the ~1e4 points become
	# a 500 dpi raster layer while axes and labels stay editable text.
	ax.scatter(
		np.log10(umi[ok]), nf[ok],
		s=3, alpha=0.3, linewidths=0, color='#106e78', rasterized=True,
		)
	ax.set_xlabel('log10(total UMI per cell)')
	ax.set_ylabel('Nuclear fraction (intronic / [intronic + exonic])')
	ax.set_title(f'{sampleid}\nn={int(ok.sum())} cells, median NF={np.nanmedian(nf[ok]):.3f}', fontsize=9)
	ax.set_ylim(0, 1)
	fig.tight_layout()

	stem = out_pdf[:-len('.pdf')]
	qcutil.savefig(fig, stem)
	plt.close(fig)


if __name__ == '__main__':
	main()
