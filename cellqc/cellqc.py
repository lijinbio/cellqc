#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

__version__ = "0.0.1"

import sys
import os
import click

def runcmd(cmd):
	print(f'Start running: {cmd}')
	try:
		exitcode=os.system(cmd)
		if exitcode!=0:
			print(f'Error: {cmd} failed.', file=sys.stderr)
			sys.exit(-1)
	except OSError as e:
		print("Execution failed: ", e, file=sys.stderr)
		sys.exit(-1)
	print(f'Finish running: {cmd}')
	return exitcode

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-c', '--configfile', required=True, type=click.Path(exists=True), help='Configuration file in YAML format.')
@click.option('-r', '--rule', type=click.STRING, help='Force to re-run a rule and its downstream. Available: soupx, dropkick, h5subset, doubletfinder.')
@click.option('-d', '--outdir', type=click.Path(), default='.', help='Outdir. (default: ".")')
@click.option('-t', '--numthreads', type=click.INT, default=4, help='Number of threads. (default: 4)')
@click.option('-D', '--dagonly', is_flag=True, help='Generate the DAG of jobs only. Useful to update dag.pdf. (default: False)')
@click.option('-R', '--reportonly', is_flag=True, help='Generate the report only. Useful to update report.html. (default: False)')
@click.option('-S', '--summaryonly', is_flag=True, help='Generate the detailed summary only. Useful to update summary.txt. (default: False)')
@click.version_option(version=__version__)
def main(configfile, rule, outdir, numthreads, dagonly, reportonly, summaryonly):
	"""
cellqc: standardized quality control pipeline of single-cell RNA-Seq data.

\b
Example:
  cellqc -c config.yaml

\b
Date: 2022/09/21
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	srcdir=os.path.dirname(os.path.abspath(__file__))
	configdir=os.path.dirname(os.path.abspath(configfile))
	cmdstr=f'snakemake -s {srcdir}/Snakefile --configfile {configfile} -C configdir={configdir} -d {outdir}'
	cmdrun=cmdstr+f' -j {numthreads}'
	cmdreport=cmdstr+f' --report report.html'
	cmddag=cmdstr+f' --dag | tee {outdir}/dag.dot | dot -Tpdf > {outdir}/dag.pdf'
	cmdsummary=cmdstr+f' -D -c1 | tee {outdir}/summary.txt'
	if rule is not None:
		cmdrun+=f' -R {rule}'

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

if __name__ == "__main__":
	main()
