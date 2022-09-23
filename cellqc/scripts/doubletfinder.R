# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DoubletFinder))
suppressPackageStartupMessages(library(ggplot2))

x=Read10X_h5(snakemake@input[[1]])
x=CreateSeuratObject(counts=x)
x=x %>% NormalizeData() %>% FindVariableFeatures(selection.method='vst', nfeatures=2000) %>% ScaleData() %>% RunPCA() %>% RunTSNE() %>% RunUMAP(dims=1:30)

outdir=snakemake@output[[1]]
if (!dir.exists(outdir)) {
	dir.create(outdir)
}

if (snakemake@params[['findpK']]) {
	sweepx=paramSweep_v3(x, PCs=1:10, sct=F, num.cores=snakemake@params[['numthreads']])
	sweepstats=summarizeSweep(sweepx, GT=F)
	pdf(sprintf('%s/findpK_bcmetric_raw.pdf', outdir), width=8, height=7.5)
	bcmetric=find.pK(sweepstats)
	dev.off()
	write.table(bcmetric, file=gzfile(sprintf('%s/findpK_bcmetric.txt.gz', outdir)), quote=F, sep='\t', row.names=F, col.names=T)
	bcmetric$pK=as.numeric(levels(bcmetric$pK))[bcmetric$pK]

	p=ggplot(data=bcmetric, aes(x=pK, y=BCmetric)) +
	geom_line(color='blue') +
	geom_point(color='blue') +
	scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
	scale_y_continuous(breaks=scales::pretty_breaks(n=10)) +
	xlab('pK') +
	ylab('Mean-variance-normalized bimodality coefficient (meanBC/varBC)') +
	ggtitle('') +
	theme_bw() +
	theme(
		plot.background=element_blank()
		, panel.grid.major=element_blank()
		, panel.grid.minor=element_blank()
		, panel.border=element_blank()
		, plot.title=element_text(hjust=0.5)
		, axis.line=element_line(color='black')
		, axis.text=element_text(color='black')
		, legend.title=element_blank()
		, legend.position='right'
		)
	ggsave(p, file=sprintf('%s/findpK_bcmetric.pdf', outdir), width=6, height=5, units='in', useDingbats=F)
	pKopt=bcmetric$pK[which.max(bcmetric[, 'BCmetric'])]
} else {
	pKopt=snakemake@params[['pK']]
}

nrun=snakemake@params[['nrun']]
doubletratio=round(ncol(x)*0.1/(nrun*13000), digits=2)
nExp_poi=round(doubletratio*ncol(x))
cat(sprintf('%s: expected %f of %d cells is %d\n', basename(snakemake@input[[1]]), doubletratio, nrow(x@meta.data), nExp_poi))

res=doubletFinder_v3(x, PCs=1:10, pN=0.25, pK=pKopt, nExp=nExp_poi, reuse.pANN=F, sct=F)

metadata=res@meta.data
header=names(metadata)
value=header[pmatch('pANN', header)]
group=header[pmatch('DF.classifications', header)]
metadata[, group]=factor(metadata[, group], c('Singlet', 'Doublet'))

p=ggplot(data=metadata, aes_q(x=as.name(group), y=as.name(value), fill=as.name(group))) +
geom_violin(trim=T, show.legend=F) +
geom_boxplot(width=0.1, fill='white', outlier.shape=NA) +
xlab('') +
ylab('pANN') +
ggtitle(value) +
theme_bw() +
theme(
	plot.background=element_blank()
	, panel.grid.major=element_blank()
	, panel.grid.minor=element_blank()
	, panel.border=element_blank()
	, plot.title=element_text(hjust=0.5)
	, axis.line=element_line(color='black')
	, axis.text=element_text(color='black')
	, legend.title=element_blank()
	, legend.position='right'
	, axis.text.x=element_text(angle=45, vjust=1, hjust=1)
	)
ggsave(p, file=sprintf('%s/pANN_violin_pK%s.pdf', outdir, pKopt), width=6, height=5, units='in', useDingbats=F)

pdf(sprintf('%s/umap_doublet_pK%s.pdf', outdir, pKopt), width=8, height=7.5)
print(DimPlot(res, reduction='umap', group.by=group))
dev.off()

pdf(sprintf('%s/tsne_doublet_pK%s.pdf', outdir, pKopt), width=8, height=7.5)
print(DimPlot(res, reduction='tsne', group.by=group))
dev.off()

res$pANN=res@meta.data[, value]
res$DF.classifications=res@meta.data[, group]
utils::write.table(
	data.frame(barcode=rownames(res@meta.data), res@meta.data)
	, file=gzfile(sprintf('%s/doubletfinder_pK%s_metadata.txt.gz', outdir, pKopt))
	, quote=F
	, sep='\t'
	, row.names=F
	, col.names=T
	)

res=subset(res, cells=rownames(res@meta.data)[res@meta.data[, group]=='Singlet'])
saveRDS(res, file=snakemake@output[[2]])
