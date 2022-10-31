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

## Filterbycounts
tmp=[pd.read_csv(f, sep='\t', header=0) for f in snakemake.input.filterbycountncell]
_filterbycountncell=pd.concat(tmp, ignore_index=True).to_html() if len(tmp)>0 else None
_filterbycountpltbf={f : data_uri_from_file(f) for f in snakemake.input.filterbycountpltbf}
_filterbycountpltaf={f : data_uri_from_file(f) for f in snakemake.input.filterbycountpltaf}

## DoubletFinder
tmp=[pd.read_csv(f, sep='\t', header=0) for f in snakemake.input.doubletratio]
_doubletratio=pd.concat(tmp, ignore_index=True).to_html() if len(tmp)>0 else None
_doubletpannviolin={f : data_uri_from_file(f) for f in snakemake.input.doubletpannviolin}
_doublettsne={f : data_uri_from_file(f) for f in snakemake.input.doublettsne}
_doubletumap={f : data_uri_from_file(f) for f in snakemake.input.doubletumap}

## scPred

## Render report file from a Jinja template
outfile=snakemake.output[0]
with open(outfile, mode="w", encoding="utf-8") as f:
	f.write(
		template.render(
			cellrangersummary=_cellrangersummary,
			soupxrhoEst=_soupxrhoEst,
			soupxrho=_soupxrho,
			filterbycountncell=_filterbycountncell,
			filterbycountpltbf=_filterbycountpltbf,
			filterbycountpltaf=_filterbycountpltaf,
			doubletratio=_doubletratio,
			doubletpannviolin=_doubletpannviolin,
			doublettsne=_doublettsne,
			doubletumap=_doubletumap,
			)
		)
