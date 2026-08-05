#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:
#
# Scenario: nreaction > 1.
#
# One Cell Ranger run can pool several 10x reactions. The expected doublet
# fraction is rate * ncell / (nreaction * capacity), so pooling n reactions
# divides the expected fraction by n -- the cells came from n separate
# emulsions, and only within-reaction collisions are doublets.
#
#   bash tests/nreaction/main.sh                             # arithmetic check only
#   CELLRANGER=/path/to/outs bash tests/nreaction/main.sh    # + full pipeline run
#
# WORKDIR defaults to $PWD; use shared storage, not /tmp, if this runs on a
# compute node (/tmp is node-local).

set -euo pipefail

WORKDIR="${WORKDIR:-$PWD}"
OUTDIR="${OUTDIR:-$WORKDIR/out_nreaction}"
NREACTION="${NREACTION:-2}"

# The claim the scenario exists to protect, checked without any data.
python - "$NREACTION" <<'PY'
import sys
from cellqc import qcutil

n = int(sys.argv[1])
ncell, rate, capacity = 13000, 0.1, 13000
one, _ = qcutil.expected_doublet_rate(ncell, 1, rate, capacity)
many, _ = qcutil.expected_doublet_rate(ncell, n, rate, capacity)
assert abs(one / n - many) < 1e-9, f'expected {one}/{n}, got {many}'
print(f'[nreaction] expected doublet fraction: {one} at nreaction=1, {many} at nreaction={n}  ok')
PY

if [[ -z ${CELLRANGER:-} ]]; then
	echo '[nreaction] CELLRANGER is not set: skipping the pipeline run.'
	echo '[nreaction] Set it to a Cell Ranger outs/ directory to run the scenario end to end.'
	exit 0
fi

cd "$WORKDIR"

printf 'sample\tcellranger\tnreaction\n' > samples.txt
printf 'pooled\t%s\t%s\n' "$CELLRANGER" "$NREACTION" >> samples.txt

cat > config.yaml <<'YAML'
seed: 42
ambient:
  method: soupx
filterbycount:
  mincount: 500
  minfeature: 300
  mito: 10
doublet:
  run: [doubletfinder, scdblfinder]
  decider: doubletfinder
  findpK: false
  pK: 0.01
  rate: 0.1
  capacity: 13000
YAML

cellqc -d "$OUTDIR" -t "${THREADS:-8}" -c config.yaml -- samples.txt

echo '[nreaction] doublet summary:'
cat "$OUTDIR"/result/*_doublet_summary.txt
