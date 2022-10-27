# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scPred))
suppressPackageStartupMessages(library(SeuratDisk))

seuratdimplotgroupby=function(data, reduct='umap', group='cluster', label=T, repel=T, outfile, width=6, height=5) {
	x=Seurat::DimPlot(data, reduction=reduct, group.by=group, label=label, repel=repel)
	if (endsWith(outfile, '.pdf')) {
		pdf(outfile, width=width, height=height)
		print(x)
		dev.off()
	} else if (endsWith(outfile, '.png')) {
		png(outfile, width=width, height=height, units='in', res=500)
		print(x)
		dev.off()
	}
}
write.txt=function(..., file=stdout(), quote=F, sep='\t', row.names=F, col.names=T) {
	if (is.character(file) && endsWith(file, '.gz')) {
		file=gzfile(file)
	}
	utils::write.table(..., file=file, quote=quote, sep=sep, row.names=row.names, col.names=col.names)
}
seuratfeatureplot=function(data, features, slot='data', color=c('blue', 'white', 'red'), outfile, width=6, height=5) {
	x=Seurat::FeaturePlot(data, features, slot=slot) + ggplot2::scale_colour_gradientn(colours=color)
	if (endsWith(outfile, '.pdf')) {
		pdf(outfile, width=width, height=height)
		print(x)
		dev.off()
	} else if (endsWith(outfile, '.png')) {
		png(outfile, width=width, height=height, units='in', res=500)
		print(x)
		dev.off()
	}
}

infile=snakemake@input[[1]]
outdir=snakemake@output[[1]]
resultfile=snakemake@output[[2]]
contingencyfile=snakemake@output[[3]]

if (!dir.exists(outdir)) {
	dir.create(outdir)
}

if (endsWith(infile, '.h5')) {
	x=Read10X_h5(infile)
	x=CreateSeuratObject(counts=x)
} else if (endsWith(infile, '.rds')) {
	x=readRDS(infile)
} else if (endsWith(infile, '.h5seurat')) {
	x=LoadH5Seurat(infile, assay='RNA')
} else {
	write('Error: wrong infile. Please input either .rds or .h5seurat file\n', stderr())
	q(status=1)
}
x=x %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunTSNE() %>% RunUMAP(dims=1:30) %>% FindNeighbors(dims=1:30) %>% FindClusters()

seuratdimplotgroupby(x, reduct='pca', group='seurat_clusters', outfile=sprintf('%s/cluster_pca.png', outdir))
seuratdimplotgroupby(x, reduct='tsne', group='seurat_clusters', outfile=sprintf('%s/cluster_tsne.png', outdir))
seuratdimplotgroupby(x, reduct='umap', group='seurat_clusters', outfile=sprintf('%s/cluster_umap.png', outdir))

ref=readRDS(snakemake@params[['ref']])
x=scPredict(x, ref, threshold=snakemake@params[['threshold']])
seuratdimplotgroupby(x, reduct='scpred', group='scpred_prediction', outfile=sprintf('%s/predict_scpred.png', outdir))
seuratdimplotgroupby(x, reduct='umap', group='scpred_prediction', outfile=sprintf('%s/predict_umap.png', outdir))

metadata=x@meta.data
write.txt(data.frame(barcode=rownames(metadata), metadata), file=sprintf('%s/metadata.txt.gz', outdir))
features=c('nCount_RNA', 'nFeature_RNA', setdiff(grep('^scpred_', names(metadata), value=T), c('scpred_prediction', 'scpred_no_rejection')))
invisible(
	lapply(features
		, function(feature) {
			seuratfeatureplot(x, features=feature, outfile=sprintf('%s/feature_%s.png', outdir, feature))
		})
	)
tmp=crossTab(x, 'seurat_clusters', 'scpred_prediction')
write.txt(cbind(cluster=rownames(tmp), tmp) , file=contingencyfile)

x=CreateSeuratObject(counts=x[['RNA']]@counts, meta.data=x@meta.data)
SaveH5Seurat(x, resultfile)
