# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import scanpy as sc
import anndata as ad
import pandas as pd

infile, metafile, sampleid, outfile=snakemake.input[0], snakemake.input[1], snakemake.params.sampleid, snakemake.output[0]

def h5ad2addmetadata(indata, metadata):
	common=indata.obs.index.intersection(metadata.index)
	if len(common)<len(indata.obs.index):
		print(f"Warning: some cells are discarded due to not in metadata {len(common)} < {len(indata.obs.index)}.")
	if len(common)==0:
		print("Warning: empty cell barcode in metadata.")
		return indata

	indata=indata[common].copy()
	metadata=metadata.loc[common]
	header=metadata.columns.tolist()
	for h in header:
		indata.obs[h]=metadata[h].tolist()
	return indata

# 0. Read .h5ad
x=sc.read_h5ad(infile)

# 1. Clean x.raw
result=ad.AnnData(X=x.X, obs=x.obs, var=x.var)
for obsm in x.obsm_keys():
	result.obsm[obsm]=x.obsm[obsm]

# 2. Add nulcear fraction
metadata=pd.read_csv(metafile, header=0, index_col=0, sep='\t')
result=h5ad2addmetadata(result, metadata)

# 3. Add prefix to cell barcode
result.obs.index=sampleid+'_'+result.obs.index

# 4. Unique var names
result.var_names_make_unique(join='.')

sc.write(filename=outfile, adata=result)
