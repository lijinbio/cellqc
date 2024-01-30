from snakemake.utils import validate
import pandas as pd
import yaml
import json
from pathlib import Path

samplefile, outdir, nowtimestr=config['samplefile'], config['outdir'], config['nowtimestr']
config['samples']=samplefile

wildcard_constraints:
	sample="[^/]+",

if config['configfile']=='None':
	config['configfile']=None

# Default parameters
default_params={
	'dropkick': {
		'skip': False,
		'method': 'multiotsu',
		'numthreads': 1,
	},
	'filterbycount': {
		'mincount': 500,
		'minfeature': 300,
		'mito': 10,
	},
	'doubletfinder': {
		'skip': False,
		'findpK': False,
		'numthreads': 5,
		'pK': 0.01,
	},
	'scpred': {
		'skip': True,
		'reference': '/path_to_reference/scPred_trainmodel_RNA_svmRadialWeights_scpred.rds',
		'threshold': 0.9,
	},
}

# Assign default params to config if not existing
for key, paramdict in default_params.items():
	if key not in config:
		config[key]=paramdict
	else:
		for subkey, value in paramdict.items():
			if subkey not in config[key]:
				config[key][subkey]=value
			else:
				print("Info: {key} and {subkey} are in config. So, skipping default assignment.")

if config['configfile']:
	configdir=Path(config['configfile']).parent
	config['scpred']['reference']=str(configdir / config['scpred']['reference'])
sampledir=str(Path(config['samples']).parent)

validate(config, schema='../schemas/config.schema.yaml')
samples=pd.read_table(config['samples']).set_index('sample', drop=False)
validate(samples, '../schemas/samples.schema.yaml')

# debug parameters
print(json.dumps(config, indent=4))
with open(f"config_{nowtimestr}.yaml", 'w') as f:
	yaml.dump(config, f, sort_keys=False)
