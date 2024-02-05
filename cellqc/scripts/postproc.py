# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import scanpy as sc
import anndata as ad

infile, sampleid, outfile=snakemake.input[0], snakemake.params.sampleid, snakemake.output[0]

x=sc.read(infile)

# 1. Clean x.raw
result=ad.AnnData(X=x.X, obs=x.obs, var=x.var)
for obsm in x.obsm_keys():
	result.obsm[obsm]=x.obsm[obsm]

# 2. Add prefix to cell barcode
result.obs.index=sampleid+'_'+result.obs.index

# 3. Unique var names
result.var_names_make_unique(join='.')

sc.write(filename=outfile, adata=result)
