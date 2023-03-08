#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug
outdir=$(mrrdir.sh)

slurmtaco.sh -n r03 -- "$(cat <<EOF
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate cellqc_v0.0.6
cellqc -n -d "$outdir" -- samples.txt
EOF
)"
