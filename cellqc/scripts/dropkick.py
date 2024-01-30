# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import sys
f, numthreads, thresh_methods=snakemake.input[0], snakemake.params.numthreads, snakemake.params.method

import scanpy as sc
sc.set_figure_params(color_map='viridis', frameon=False)

import os
if os.path.isdir(f):
	x=sc.read_10x_mtx(f)
elif f.endswith('.h5'):
	x=sc.read_10x_h5(f)
elif f.endswith('.h5ad'):
	x=sc.read(f)
x.var_names_make_unique(join='.')
print(vars(x))

import dropkick as dk
x=dk.recipe_dropkick(x, n_hvgs=None, X_final='raw_counts')
dk.qc_summary(x, save_to=snakemake.output[0])
print(vars(x))

x_classifier=dk.dropkick(x, n_jobs=numthreads, thresh_methods=thresh_methods)

dk.score_plot(x, save_to=snakemake.output[1])
dk.coef_inventory(x)
dk.coef_plot(x, save_to=snakemake.output[2])

x.obs.dropkick_label=x.obs.dropkick_label.astype('string') # A potential bug without this conversion
sc.write(filename=snakemake.output[3], adata=x)
x.obs['barcode']=x.obs.index
x.obs.to_csv(snakemake.output[4], sep='\t', index=False)
x.var['symbol']=x.var.index
x.var.to_csv(snakemake.output[5], sep='\t', index=False)
