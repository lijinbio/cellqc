# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(SeuratDisk))

infile=snakemake@input[[1]]
outdir=snakemake@output[[1]]
mincount=snakemake@params[['mincount']]
minfeature=snakemake@params[['minfeature']]
mito=snakemake@params[['mito']]
sampleid=snakemake@params[['sampleid']]
outfile=snakemake@output[[2]]
statfile=snakemake@output[[3]]

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

if (!dir.exists(outdir)) {
	dir.create(outdir)
}

x[['percent.mt']]=PercentageFeatureSet(x, pattern='^MT-|^mt-')

Idents(x)=sampleid
p=VlnPlot(x, features=c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol=3, log=F, pt.size=0) + NoLegend()
ggsave(p, file=sprintf('%s/feature_bf.pdf', outdir), width=9, height=7.5, units='in', useDingbats=F)

ncell_bf=ncol(x)
x=subset(x, subset=nCount_RNA>=mincount & nFeature_RNA>=minfeature & percent.mt<=mito)
ncell_af=ncol(x)

p=VlnPlot(x, features=c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol=3, log=F, pt.size=0) + NoLegend()
ggsave(p, file=sprintf('%s/feature_af.pdf', outdir), width=9, height=7.5, units='in', useDingbats=F)

utils::write.table(data.frame(sampleid=sampleid, ncell_bf=ncell_bf, ncell_af=ncell_af), file=statfile, quote=F, sep='\t', row.names=F, col.names=T)
SaveH5Seurat(x, outfile)
