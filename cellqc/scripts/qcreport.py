# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import os
import sys
import pandas as pd
from pathlib import Path
import mimetypes
import base64
from jinja2 import Environment, FileSystemLoader

# From snakemake/report/__init__.py
def data_uri(data, filename, encoding="utf8", mime="text/plain"):
	"""Craft a base64 data URI from file with proper encoding and mimetype."""
	data=base64.b64encode(data)
	uri="data:{mime};charset={charset};filename={filename};base64,{data}" "".format(
		filename=filename, mime=mime, charset=encoding, data=data.decode("utf-8")
		)
	return uri

def mime_from_file(file):
	mime, encoding=mimetypes.guess_type(file)
	if mime is None:
		mime="text/plain"
		print(f"Could not detect mimetype for {file}, assuming text/plain.", file=sys.stderr)
	return mime, encoding

def data_uri_from_file(file, defaultenc="utf8"):
	"""Craft a base64 data URI from file with proper encoding and mimetype."""
	if isinstance(file, Path):
		file=str(file)
	mime, encoding=mime_from_file(file)
	if encoding is None:
		encoding=defaultenc
	with open(file, "rb") as f:
		return data_uri(f.read(), os.path.basename(file), encoding, mime)

# From snakemake/report/data/common.py
def get_resource_as_string(path):
	return open(Path(__file__).parent / "template" / path).read()

env=Environment(loader=FileSystemLoader(Path(__file__).parent / "template"))
env.filters["get_resource_as_string"]=get_resource_as_string
template=env.get_template("index.html.jinja2")

## Cellranger metrics summary
tmp=[]
for k, v in snakemake.params.samples['cellranger'].to_dict().items():
	vf=os.path.join(snakemake.params.sampledir, f"{v}/metrics_summary.csv")
	if os.path.exists(vf):
		x=pd.read_csv(vf)
		x.insert(0, 'sampleid', k)
		tmp+=[x]
_cellrangersummary=pd.concat(tmp, ignore_index=True).to_html() if len(tmp)>0 else None

## SoupX
tmp=[pd.read_csv(f, sep='\t', header=0) for f in snakemake.input.soupxrhoEst]
_soupxrhoEst=pd.concat(tmp, ignore_index=True).to_html() if len(tmp)>0 else None
_soupxrho={f : data_uri_from_file(f) for f in snakemake.input.soupxrho}

## Dropkick
tmp=[pd.read_csv(f, sep='\t', header=0) for f in snakemake.input.dropkickstat]
_dropkickstat=pd.concat(tmp, ignore_index=True).to_html() if len(tmp)>0 else None
_dropkicksummary={f : data_uri_from_file(f) for f in snakemake.input.dropkicksummary} if len(snakemake.input.dropkicksummary)>0 else None
_dropkickscore={f : data_uri_from_file(f) for f in snakemake.input.dropkickscore} if len(snakemake.input.dropkickscore)>0 else None

## Filterbycounts
tmp=[pd.read_csv(f'{dir}/filter_ncell.txt', sep='\t', header=0) for dir in snakemake.input.filterbycountdir]
_filterbycountncell=pd.concat(tmp, ignore_index=True).to_html() if len(tmp)>0 else None
_filterbycountpltbf={dir : data_uri_from_file(f'{dir}/feature_bf.pdf') for dir in snakemake.input.filterbycountdir}
_filterbycountpltaf={dir : data_uri_from_file(f'{dir}/feature_af.pdf') for dir in snakemake.input.filterbycountdir}

## DoubletFinder
pK=snakemake.params.pK
tmp=[pd.read_csv(f'{dir}/doublet_ratio.txt', sep='\t', header=0) for dir in snakemake.input.doubletfinderdir]
_doubletratio=pd.concat(tmp, ignore_index=True).to_html() if len(tmp)>0 else None
_doubletpannviolin={dir : data_uri_from_file(f'{dir}/pANN_violin_pK{pK}.pdf') for dir in snakemake.input.doubletfinderdir}
_doublettsne={dir : data_uri_from_file(f'{dir}/tsne_doublet_pK{pK}.pdf') for dir in snakemake.input.doubletfinderdir}
_doubletumap={dir : data_uri_from_file(f'{dir}/umap_doublet_pK{pK}.pdf') for dir in snakemake.input.doubletfinderdir}

## scPred
_scpredcontingency={dir : pd.read_csv(f'{dir}/contingency.txt', sep='\t', header=0).to_html() for dir in snakemake.input.scpreddir} if len(snakemake.input.scpreddir)>0 else None
_scpredfeaturemax={dir : data_uri_from_file(f'{dir}/feature_scpred_max.png') for dir in snakemake.input.scpreddir} if len(snakemake.input.scpreddir)>0 else None
_scpredpredictumap={dir : data_uri_from_file(f'{dir}/predict_umap.png') for dir in snakemake.input.scpreddir} if len(snakemake.input.scpreddir)>0 else None

## Render report file from a Jinja template
outfile=snakemake.output[0]
with open(outfile, mode="w", encoding="utf-8") as f:
	f.write(
		template.render(
			cellrangersummary=_cellrangersummary,
			soupxrhoEst=_soupxrhoEst,
			soupxrho=_soupxrho,
			dropkickstat=_dropkickstat,
			dropkicksummary=_dropkicksummary,
			dropkickscore=_dropkickscore,
			filterbycountncell=_filterbycountncell,
			filterbycountpltbf=_filterbycountpltbf,
			filterbycountpltaf=_filterbycountpltaf,
			doubletratio=_doubletratio,
			doubletpannviolin=_doubletpannviolin,
			doublettsne=_doublettsne,
			doubletumap=_doubletumap,
			scpredcontingency=_scpredcontingency,
			scpredfeaturemax=_scpredfeaturemax,
			scpredpredictumap=_scpredpredictumap,
			)
		)
