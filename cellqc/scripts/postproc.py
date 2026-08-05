# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:
"""Post-processing for downstream integration. Writes the final matrix.

Prefixes the cell barcode with the sample ID, makes var names unique, drops the
`raw` layer, and attaches the nuclear fraction when it was computed. The nuclear
fraction input is absent for samples with no Cell Ranger BAM (detected at DAG
construction), so its presence is checked rather than assumed.

This is the last stage, so its output is `result/{sample}.h5ad` -- the matrix a
user takes away -- alongside `.obs` and `.var` as gzipped TSVs.
"""

import anndata as ad
import pandas as pd
import scanpy as sc

from cellqc import qcutil

infile = snakemake.input['h5ad']
metafile = snakemake.input.get('nf', None)
sampleid = snakemake.params['sampleid']
outfile = snakemake.output['h5ad']
out_obs = snakemake.output['obs']
out_var = snakemake.output['var']


def add_metadata(indata, metadata, what):
	common = indata.obs.index.intersection(metadata.index)
	if len(common) == 0:
		raise ValueError(
			f'{sampleid}: no barcode overlap between the matrix and {what}. '
			'Refusing to write a matrix with a silently empty annotation.'
			)
	if len(common) < len(indata.obs.index):
		# Report rather than drop: subsetting the matrix to whatever an
		# annotation happens to cover is how cells disappear unnoticed.
		print(
			f'[postproc] {sampleid}: WARNING {what} covers {len(common)}/{len(indata.obs.index)} '
			'cells; the remainder get NaN and are KEPT',
			flush=True,
			)
	metadata = metadata.reindex(indata.obs.index)
	for h in metadata.columns:
		indata.obs[h] = metadata[h].to_numpy()
	return indata


def main():
	x = sc.read_h5ad(infile)
	n_in = x.n_obs

	# Clean x.raw -- carrying a stale raw layer through integration is a common
	# source of silently mismatched matrices.
	result = ad.AnnData(X=x.X, obs=x.obs.copy(), var=x.var.copy(), uns=dict(x.uns))
	for obsm in x.obsm_keys():
		result.obsm[obsm] = x.obsm[obsm]

	if metafile:
		metadata = pd.read_csv(metafile, header=0, index_col=0, sep='\t')
		result = add_metadata(result, metadata, 'the nuclear fraction table')
		result.uns['cellqc_nuclear_fraction'] = True
	else:
		# Recorded explicitly so a missing column is never mistaken downstream
		# for a computed-and-zero one.
		result.uns['cellqc_nuclear_fraction'] = False
		print(
			f'[postproc] {sampleid}: no nuclear fraction (no Cell Ranger BAM for this sample)',
			flush=True,
			)

	result.obs.index = sampleid + '_' + result.obs.index.astype(str)
	result.obs['sampleid'] = sampleid
	result.var_names_make_unique(join='.')

	assert result.n_obs == n_in, f'{sampleid}: postproc changed the cell count {n_in} -> {result.n_obs}'
	result.write_h5ad(outfile, compression='gzip')
	qcutil.write_obs_var(result, out_obs, out_var)
	print(f'[postproc] {sampleid}: {result.n_obs} cells x {result.n_vars} genes written', flush=True)


if __name__ == '__main__':
	main()
