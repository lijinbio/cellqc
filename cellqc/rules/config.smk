from snakemake.utils import validate
import pandas as pd
import yaml
import json
import os
import sys
from pathlib import Path

# NB: do not `import cellqc` here. Snakemake puts the Snakefile's directory on
# sys.path, where cellqc/cellqc.py shadows the installed `cellqc` package, so the
# import resolves to the CLI module and fails. The `script:` steps run in their
# own process and can import the package normally; .smk files keep their own
# copies of the few helpers they need.
BAM_NAME='possorted_genome_bam.bam'


def has_bam(cellranger_dir):
	bam=os.path.join(cellranger_dir, BAM_NAME)
	return os.path.exists(bam) and os.path.exists(bam+'.bai')

samplefile, outdir, nowtimestr=config['samplefile'], config['outdir'], config['nowtimestr']
config['samples']=samplefile

wildcard_constraints:
	sample="[^/]+",

if config['configfile']=='None':
	config['configfile']=None

# Default parameters
default_params={
	'seed': 42,
	'ambient': {
		'method': 'soupx',
		'compare': [],
	},
	'nuclear_fraction': {
		'numthreads': 12,
		'cbtag': 'CB',
		'retag': 'RE',
		'exontag': 'E',
		'introntag': 'N',
	},
	'filterbycount': {
		'mincount': 500,
		'minfeature': 300,
		'mito': 10,
	},
	# QC gene sets. Kept here rather than in cellqc/qcutil.py's GENE_SETS so the
	# effective definition lands in config_<timestamp>.yaml with everything else:
	# which genes counted as mitochondrial is not reconstructable from the matrix
	# afterwards, and a run whose reference names them differently is exactly the
	# run whose numbers a reader will question. The two must stay in step -- see
	# the note in qcutil.py for what `patterns`, `symbols` and `exclude` mean.
	# Only 'mt' is filtered on (filterbycount.mito); the rest are recorded.
	'geneset': {
		'mt': {
			'label': '% mitochondrial',
			'patterns': [r'^MT-'],
			'symbols': [
				'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6',
				'COX1', 'COX2', 'COX3', 'ATP6', 'ATP8', 'CYTB',
			],
			'exclude': [],
		},
		'ribo': {
			'label': '% ribosomal',
			'patterns': [r'^RP[SL]\d', r'^RPLP\d', r'^RPSA$'],
			'symbols': [],
			'exclude': [r'^RPS6K', r'^RPS19BP'],
		},
		'hb': {
			'label': '% hemoglobin',
			'patterns': [r'^HB[ABDEGMQZ]([0-9][AB]?)?$', r'^HB[AB]-[A-Z0-9]+$'],
			'symbols': [],
			'exclude': [],
		},
	},
	'doublet': {
		'run': ['doubletfinder', 'scdblfinder'],
		'decider': 'doubletfinder',
		'findpK': False,
		'numthreads': 5,
		'pK': 0.01,
		'rate': 0.1,
		'capacity': 13000,
		'nreaction': 1,
	},
}

# v0.1.0 sections that no longer exist. Warn and drop rather than fail, so old
# config files keep running instead of forcing an edit before a re-run.
removed_sections={
	'dropkick': 'empty-droplet calling was removed in v0.2.0; Cell Ranger EmptyDrops is used alone',
	'scpred': 'cell-type annotation was removed in v0.2.0; annotate downstream of cellqc',
}
for key, why in removed_sections.items():
	if key in config:
		print(f"Warning: config section '{key}' is obsolete and will be ignored ({why}).", file=sys.stderr)
		del config[key]

# v0.1.0 called this section 'doubletfinder'; it is 'doublet' now because it
# covers more than one caller. Migrate the shared keys so old configs still work.
if 'doubletfinder' in config:
	print("Warning: config section 'doubletfinder' was renamed to 'doublet' in v0.2.0; migrating.", file=sys.stderr)
	legacy=config.pop('doubletfinder')
	migrated=config.setdefault('doublet', {})
	for k, v in legacy.items():
		migrated.setdefault(k, v)

# There is no doublet skip flag, for the same reason the nuclear-fraction step
# has none: which callers run is `doublet.run`, and a caller that should not run
# is left out of it. `skip: false` was redundant with that; `skip: true` has no
# equivalent, so it is an error rather than a warning -- silently ignoring it
# would write a matrix with doublets removed to a config that asked for none.
if 'skip' in config.get('doublet', {}):
	if config['doublet'].pop('skip'):
		raise ValueError(
			"config key 'doublet.skip' was removed: doublet detection always runs. "
			"There is no configuration that keeps every called doublet, so a run that "
			"previously set 'doublet.skip: true' would now silently lose cells. "
			"Remove the key and choose the callers with doublet.run/doublet.decider."
			)
	print(
		"Warning: config key 'doublet.skip' is obsolete and was ignored; select callers with doublet.run.",
		file=sys.stderr)

# Assign default params to config if not existing
for key, paramdict in default_params.items():
	if not isinstance(paramdict, dict):
		config.setdefault(key, paramdict)
		continue
	if key not in config:
		config[key]=paramdict
	else:
		for subkey, value in paramdict.items():
			if subkey not in config[key]:
				config[key][subkey]=value

# The generic merge above is one level deep, which is the rule for `geneset`
# too: naming a set in the config replaces that set's definition outright rather
# than merging key by key. A half-merged set is the worse behaviour -- someone
# who writes their own `patterns` for a reference and inherits the default
# `symbols` gets the macaque fallback silently applied to it.
if config['geneset'].get('mt') is None:
	raise ValueError(
		"config section 'geneset' must define 'mt': filterbycount.mito is a threshold on "
		"pct_counts_mt, so the mitochondrial gene set cannot be removed."
		)

# The deciding caller must actually run, otherwise its metadata never exists.
if config['doublet']['decider'] not in config['doublet']['run']:
	raise ValueError(
		f"doublet.decider={config['doublet']['decider']!r} is not in doublet.run={config['doublet']['run']!r}. "
		"The caller that removes cells must be one of the callers that runs."
		)

if config['ambient']['method'] in config['ambient']['compare']:
	print(
		f"Warning: ambient.method={config['ambient']['method']!r} also listed in ambient.compare; "
		"dropping the duplicate.", file=sys.stderr)
	config['ambient']['compare']=[m for m in config['ambient']['compare'] if m!=config['ambient']['method']]

sampledir=str(Path(config['samples']).parent)

validate(config, schema='../schemas/config.schema.yaml')
samples=pd.read_table(config['samples']).set_index('sample', drop=False)

# `nreaction` is per sample because a cohort can mix chemistries, but usually
# does not, so doublet.nreaction supplies the value for every sample that does
# not name its own. Filled before validation so a partly-filled column is a
# defaulted column rather than a type error.
if 'nreaction' not in samples.columns:
	samples['nreaction']=config['doublet']['nreaction']
else:
	samples['nreaction']=samples['nreaction'].fillna(config['doublet']['nreaction']).astype(int)

validate(samples, '../schemas/samples.schema.yaml')

if samples.index.duplicated().any():
	dup=sorted(set(samples.index[samples.index.duplicated()]))
	raise ValueError(f"Duplicate sample IDs in {config['samples']}: {dup}")

# The nuclear fraction needs the Cell Ranger BAM. Rather than a skip flag, detect
# it per sample at DAG-construction time: present -> the step runs, absent -> it
# is dropped for that sample. Mixed cohorts are therefore fine, and which samples
# have it is recorded so a missing column is never mistaken for a computed one.
samples['cellrangerdir']=[str(Path(sampledir) / p) for p in samples['cellranger']]
samples['has_bam']=[has_bam(p) for p in samples['cellrangerdir']]

nf_samples=samples.index[samples['has_bam']].tolist()
nobam=samples.index[~samples['has_bam']].tolist()
if nobam:
	print(
		f"Info: no indexed {'possorted_genome_bam.bam'} for {len(nobam)}/{len(samples)} sample(s): {nobam}. "
		"Skipping the nuclear-fraction step for these; all other steps run normally.",
		file=sys.stderr,
		)

missing=[s for s in samples.index if not Path(samples.loc[s, 'cellrangerdir']).is_dir()]
if missing:
	raise FileNotFoundError(
		f"Cell Ranger directory not found for sample(s) {missing}. "
		f"Paths in {config['samples']} are resolved relative to {sampledir}."
		)

# debug parameters
print(json.dumps(config, indent=4))
with open(f"config_{nowtimestr}.yaml", 'w') as f:
	yaml.dump(config, f, sort_keys=False)
