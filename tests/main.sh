#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug
outdir=$(mrrdir.sh)
slurmtaco.sh -t 4 -m 20G -- cellqc -c config.yaml -d "$outdir"
