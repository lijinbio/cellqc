#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

__version__ = "0.0.1"

import sys
import os
def run(args):
	configfile, numthreads, outdir=args.configfile.name, args.numthreads, args.outdir
	srcdir=os.path.dirname(os.path.abspath(__file__))
	configdir=os.path.dirname(os.path.abspath(configfile))
	cmd=f'snakemake -j {numthreads} -s {srcdir}/Snakefile --configfile {configfile} -C configdir={configdir} -d {outdir} --report report.html'
	if args.rule is not None:
		cmd+=f' -R {args.rule}'
	os.system(cmd)

import argparse
def main():
	parser = argparse.ArgumentParser(
		description='cellqc: standardized quality control pipeline of single-cell RNA-Seq data'
		, formatter_class=argparse.RawDescriptionHelpFormatter
		, epilog='''
Example:
  cellqc -c config.yaml

Date: 2021/03/11
Authors: Jin Li <lijin.abc@gmail.com>
'''
		)
	parser.add_argument('-c', '--configfile', type=argparse.FileType('r'), required=True, help='Configuration file in YAML format.')
	parser.add_argument('-r', '--rule', type=str, help='Force to re-run a rule and its downstream. Available: soupx, dropkick, h5subset, doubletfinder.')
	parser.add_argument('-d', '--outdir', type=str, default='.', help='Outdir. (default: ".")')
	parser.add_argument('-t', '--numthreads', type=int, default=4, help='Number of threads. (default: 4)')
	parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
	args = parser.parse_args()
	run(args)

if __name__ == "__main__":
	main()
