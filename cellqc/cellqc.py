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

\b
```
sample<TAB>cellranger[<TAB>nreaction]
```

\b
Example:
  # 1. To generate an example config.yaml by -n|--dryrun
  outdir=$(mktemp -d -u)
  cellqc -d "$outdir" -t 8 -n -- samples.txt # save default parameters to outdir/config.yaml
  ## copy from outdir/config.yaml and edit config.yaml
  cellqc -d "$outdir" -t 8 -c config.yaml -- samples.txt

\b
  # 2. To overwrite samplefile by -D|--define
  define=(
  sample:=:BCM_22_0769_RNA_fovea
  cellranger:=:/storage/chenlab/Users/jinli/wkfl/atlashumanprj/application/HCA/Retina/preproc/fixversion/v7.0.1/cellrangercount/BCM_22_0769_RNA_fovea/BCM_22_0769_RNA_fovea/outs
  nreaction:=:1
  )
  cellqc -d "$outdir" -t 8 -c config.yaml $(basharr2cmdopts.sh -o -D -- "${define[@]}") 

\b
Note:
  1. Please refer to the usage in the short tutorial available at https://github.com/lijinbio/cellqc.
  2. The "-D|--define" option can be utilized to define an individual sample by overwriting the SAMPLEFILE.

\b
Date: 2024/01/30
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
		f"-s {absdir}/Snakefile",
		f"-d {outdir}",
		f"-j {numthreads}",
		f"-C samplefile='{samplefile}' outdir='{outdir}' configfile='{configfile}' nowtimestr='{nowtimestr}'",
		]

	if configfile:
		cmdstr+=[f"--configfile {configfile}"]

	if dryrun:
		cmdstr+=[f"-n -p"]
		runcmd(cmdstr)

	else:
		cmdrun=cmdstr+[
			f"-r -p --debug-dag",
			f"--stats Snakefile_{nowtimestr}.stats",
			]
		runcmd(cmdrun)

	return 0

if __name__ == "__main__":
	sys.exit(main())
