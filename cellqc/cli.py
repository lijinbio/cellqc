#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import sys
import click
import datetime
import shutil
import subprocess
import os
import re
import pandas as pd
from pathlib import Path

def runcmd(cmdstr):
	try:
		cmdstr=' '.join(cmdstr)

		click.echo(f'Start running: {cmdstr}', file=sys.stderr)
		exitcode=subprocess.call(cmdstr, shell=True, executable=shutil.which('bash'))
		click.echo(f'Finish running: {cmdstr}', file=sys.stderr)

		if exitcode!=0:
			click.echo(f'Error: {cmdstr} failed. Exit code: {exitcode}', file=sys.stderr)
			sys.exit(-1)

	except OSError as e:
		click.echo("Execution failed: ", e, file=sys.stderr)
		sys.exit(-1)

	return exitcode

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(resolve_path=True), default='.', show_default=True, help='Outdir.')
@click.option('-c', '--configfile', type=click.Path(exists=False, resolve_path=True), help='Configuration file in the YAML format.')
@click.option('-t', '--numthreads', type=click.INT, default=4, show_default=True, help='Number of threads.')
@click.option('-n', '--dryrun', is_flag=True, help='Dry-run.')
@click.option('-D', '--define', type=click.STRING, multiple=True, help='Define KEY=VALUE pairs.')
@click.argument('samplefile', type=click.Path(exists=False, resolve_path=True), required=False)
@click.version_option()
def main(configfile, outdir, numthreads, dryrun, define, samplefile):
	"""
cellqc: a quality control pipeline of single-cell RNA-Seq data.

SAMPLEFILE is a tab-delimited file with headers, providing information for samples in the following format.
Relative "cellranger" paths are resolved against the SAMPLEFILE's own directory.

\b
```
sample<TAB>cellranger[<TAB>nreaction]
S1<TAB>/path/to/cellranger/S1/outs<TAB>1
```

\b
Example:
  # 1. A single sample, no sample file needed (-D writes one for you)
  cellqc -d out -t 8 -D sample:=:S1 -D cellranger:=:/path/to/cellranger/S1/outs

\b
  # 2. A cohort, from a sample file
  cellqc -d out -t 8 -- samples.txt                 # default parameters
  cellqc -d out -t 8 -c config.yaml -- samples.txt  # customized parameters

\b
  # 3. Dry run first: prints the jobs and writes the fully resolved
  #    configuration, defaults included, to outdir/config_<timestamp>.yaml.
  #    Copy that file, edit it, and pass it back with -c.
  cellqc -d out -t 8 -n -- samples.txt

\b
Output (per sample, under -d|--outdir):
  result/{sample}.h5ad          the final QC'd matrix, ready for integration
  result/{sample}_obs.txt.gz    .obs as a TSV, indexed by barcode
  result/{sample}_var.txt.gz    .var as a TSV, indexed by gene
  result/report.html            self-contained QC report
  result/report_slides.pdf      presentation-ready slide deck

\b
Note:
  1. Please refer to the usage in the short tutorial available at https://github.com/lijinbio/cellqc.
  2. "-D|--define" defines a single sample instead of a SAMPLEFILE. It writes
     outdir/samples_<timestamp>.txt, so the "cellranger" path must be absolute
     (a relative one would resolve against the outdir). Give it once per field:
     -D sample:=:S1 -D cellranger:=:/abs/path/outs -D nreaction:=:1
  3. Ambient RNA and doublet methods are chosen in the configuration file, not
     on the command line; see the README.

\b
Date: 2026/08/05
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	nowtimestr=datetime.datetime.now().strftime('%y%m%d_%H%M%S')
	Path(outdir).mkdir(parents=True, exist_ok=True)
	absdir=Path(__file__).parent
	os.chdir(outdir)

	if define: # overwrite samplefile
		sampleconfig={}
		for kv in define:
			m=re.match(r"^(\w+):=:(\S+)$", kv)
			if m:
				key, value=m.groups()
				sampleconfig[key]=value
			else:
				print(f"Warning: {kv} not recognized.")

		if len(sampleconfig)>0: # overwrite the sample file
			samplefile=f"{outdir}/samples_{nowtimestr}.txt"
			pd.DataFrame([sampleconfig]).to_csv(samplefile, sep='\t', index=False)

	cmdstr=[
		f"snakemake",
		f"--snakefile {absdir}/Snakefile",
		f"--directory {outdir}",
		f"--cores {numthreads}",
		f"--jobs {numthreads}",
		f"--config samplefile='{samplefile}' outdir='{outdir}' configfile='{configfile}' nowtimestr='{nowtimestr}'",
		]

	if configfile:
		cmdstr+=[f"--configfile {configfile}"]

	if dryrun:
		cmdstr+=[f"--dry-run --printshellcmds"]
		runcmd(cmdstr)

	else:
		# --use-conda is gone: v0.2.0 runs from a single environment
		# (envs/cellqc.yaml), so there are no per-rule conda directives for it to
		# act on and it only slowed startup.
		cmdstr+=[
			f"--printshellcmds --skip-script-cleanup",
			]
		runcmd(cmdstr)

	return 0

if __name__ == "__main__":
	sys.exit(main())
