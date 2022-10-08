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

cellrangersummary, outfile=snakemake.input.cellrangersummary, snakemake.output[0]
tmp=[]
for k, v in cellrangersummary:
	x=pd.read_csv(v)
	x.insert(0, 'sampleid', k)
	tmp+=[x]
cellrangersummary=pd.concat(tmp, ignore_index=True)

with open(outfile, mode="w", encoding="utf-8") as f:
	f.write(
		template.render(
			cellrangersummary=cellrangersummary.to_html(),
			)
		)
