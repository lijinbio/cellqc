#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

__version__ = "0.0.1"

import os
import click
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-c', '--configfile', required=True, type=click.Path(exists=True), help='Configuration file in YAML format.')
@click.option('-r', '--rule', type=click.STRING, help='Force to re-run a rule and its downstream. Available: soupx, dropkick, h5subset, doubletfinder.')
@click.option('-d', '--outdir', type=click.Path(), default='.', help='Outdir. (default: ".")')
@click.option('-t', '--numthreads', type=click.INT, default=4, help='Number of threads. (default: 4)')
@click.version_option(version=__version__)
def main(configfile, rule, outdir, numthreads):
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
	if rule is not None:
		cmdrun+=f' -R {rule}'
	os.system(cmdrun)
	os.system(cmdreport)

if __name__ == "__main__":
	main()
