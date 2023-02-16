from snakemake.utils import validate
import pandas as pd
import yaml
import os
import sys
from pathlib import Path

samplefile, outdir, nowtimestr=config['samplefile'], config['outdir'], config['nowtimestr']
config['samples']=samplefile

# Default parameters
if 'dropkick' not in config:
	config['dropkick']={
		'skip': False,
		'method': 'multiotsu',
		'numthreads': 1,
	}

if 'filterbycount' not in config:
	config['filterbycount']={
		'mincount': 500,
		'minfeature': 300,
		'mito': 10,
	}

if 'doubletfinder' not in config:
	config['doubletfinder']={
		'findpK': False,
		'numthreads': 5,
		'pK': 0.01,
	}

if 'scpred' not in config:
	config['scpred']={
		'skip': True,
		'reference': '/path_to_reference/scPred_trainmodel_RNA_svmRadialWeights_scpred.rds',
		'threshold': 0.9,
	}

if config['configfile']!='None':
	configdir=Path(config['configfile']).parent
	config['scpred']['reference']=str(configdir / config['scpred']['reference'])
sampledir=str(Path(config['samples']).parent)

validate(config, schema='schemas/config.schema.yaml')
samples=pd.read_table(config['samples']).set_index('sample', drop=False)
validate(samples, 'schemas/samples.schema.yaml')

# debug parameters
nowtimestr=config['nowtimestr']
with open(f"config_{nowtimestr}.yaml", 'w') as f:
	yaml.dump(config, f, sort_keys=False)
