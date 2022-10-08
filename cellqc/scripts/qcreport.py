# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

from jinja2 import Environment, FileSystemLoader
import pandas as pd
import os
# From snakemake/report/data/common.py
from pathlib import Path

def get_resource_as_string(path):
	return open(Path(__file__).parent / "template" / path).read()

env=Environment(loader=FileSystemLoader(Path(__file__).parent / "template"))
env.filters["get_resource_as_string"] = get_resource_as_string
template=env.get_template("index.html.jinja2")

## Cellranger metrics summary
tmp=[]
for k, v in snakemake.params.samples['cellranger'].to_dict().items():
	vf=os.path.join(snakemake.params.sampledir, f"{v}/metrics_summary.csv")
	if os.path.exists(vf):
		x=pd.read_csv(vf)
		x.insert(0, 'sampleid', k)
		tmp+=[x]
cellrangersummary=pd.concat(tmp, ignore_index=True).to_html() if len(tmp)>0 else None

## SoupX

## Filterbycounts

## DoubletFinder

## Render report file from a Jinja template
outfile=snakemake.output[0]
with open(outfile, mode="w", encoding="utf-8") as f:
	f.write(
		template.render(
			cellrangersummary=cellrangersummary,
			)
		)
