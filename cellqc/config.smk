from snakemake.utils import validate, min_version
min_version("7.0")
import pandas as pd
import yaml
import os
import sys
from pathlib import Path

# Default parameters
if 'samples' not in config:
	config['samples']='samples.txt'

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

if 'configfile' in config:
	configdir=Path(config['configfile']).parent
	config['samples']=str(configdir / config['samples'])
	config['scpred']['reference']=str(configdir / config['scpred']['reference'])
validate(config, schema='schemas/config.schema.yaml')

samples=pd.read_table(config['samples']).set_index('sample', drop=False)
validate(samples, 'schemas/samples.schema.yaml')

if not Path(config['samples']).exists():
	print(f"Error: {config['samples']} is not found.")
	sys.exit(-1)
sampledir=str(Path(config['samples']).parent)

# debug parameters
nowtimestr=config['nowtimestr']
configfile=f"config_{nowtimestr}.yaml"
with open(configfile, 'w') as f:
	yaml.dump(config, f, sort_keys=False)
