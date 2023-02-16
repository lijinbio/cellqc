#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import sys
import click
import datetime
import shutil
import subprocess
import os
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
@click.option('-d', '--outdir', type=click.Path(), default='.', help='Outdir. (default: ".")')
@click.option('-c', '--configfile', type=click.Path(exists=False, resolve_path=True), help='Configuration file in the YAML format.')
@click.option('-r', '--rule', type=click.STRING, help='Force to re-run a rule and its downstream. Available: soupx, dropkick, h5subset, doubletfinder.')
@click.option('-t', '--numthreads', type=click.INT, default=4, help='Number of threads. (default: 4)')
@click.option('-D', '--dagonly', is_flag=True, help='Generate the DAG of jobs only. Useful to update dag.pdf.')
@click.option('-R', '--reportonly', is_flag=True, help='Generate the report only. Useful to update report.html.')
@click.option('-S', '--summaryonly', is_flag=True, help='Generate the detailed summary only. Useful to update summary.txt.')
@click.option('-n', '--dryrun', is_flag=True, help='Dry-run.')
@click.argument('samplefile', type=click.Path(exists=False, resolve_path=True))
@click.version_option()
def main(configfile, rule, outdir, numthreads, dagonly, reportonly, summaryonly, dryrun, samplefile):
	"""
cellqc: standardized quality control pipeline of single-cell RNA-Seq data.

SAMPLEFILE is a headed tab-delimited sample file for samples in the below format.

\b
```
sample<TAB>cellranger[<TAB>nreaction]
```

\b
Example:
  outdir=$(mktemp -d -u)
  cellqc -d "$outdir" -t 8 -n -- samples.txt # save default parameters to outdir/config.yaml
  ## copy from outdir/config.yaml and edit config.yaml
  cellqc -d "$outdir" -t 8 -c config.yaml -- samples.txt

\b
Note:
  1. The directory of the `cellranger` in SAMPLEFILE should exist.

\b
Date: 2023/02/15
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	nowtimestr=datetime.datetime.now().strftime('%y%m%d_%H%M%S')
	Path(outdir).mkdir(parents=True, exist_ok=True)
	absdir=Path(__file__).parent
	os.chdir(outdir)

	cmdstr=[
		f"snakemake",
		f"-s {absdir}/Snakefile",
		f"-d {outdir}",
		f"-j {numthreads}",
		f"-C samplefile='{samplefile}' outdir='{outdir}' configfile='{configfile}' nowtimestr='{nowtimestr}'",
		]
	if rule:
		cmdstr+=[f"-R {rule}"]
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
		cmdreport=cmdstr+[
			f"--report report_{nowtimestr}.html",
			]
		cmddag=cmdstr+[
			f"--dag | tee dag_{nowtimestr}.dot | dot -Tpdf > dag_{nowtimestr}.pdf",
			]
		cmdsummary=cmdstr+[
			f"-D -c1 | tee summary_{nowtimestr}.txt",
			]

		if dagonly:
			runcmd(cmddag)
		elif reportonly:
			runcmd(cmdreport)
		elif summaryonly:
			runcmd(cmdsummary)
		else:
			runcmd(cmdrun)
			runcmd(cmddag)
			runcmd(cmdreport)
			runcmd(cmdsummary)

	return 0

if __name__ == "__main__":
	sys.exit(main())
