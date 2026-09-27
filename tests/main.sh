#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:
#
# End-to-end test of per-sample failure tolerance on a real cohort.
#
#   bash tests/main.sh [-n] SAMPLESHEET OUTDIR [THREADS]
#
#   -n           dry run: build the DAG and print the jobs
#   SAMPLESHEET  tab-delimited, with columns `sample` and `cellranger`
#                (`sampleid` and `cellranger_dir` are accepted too)
#   OUTDIR       output directory; re-running the same command resumes
#   THREADS      default 32
#
# Use a cohort that includes low-quality libraries: ones where SoupX finds no
# marker genes, or where filtering leaves too few cells for doublet detection.
# Pass criteria (see tests/README.md): the run exits 0; failing samples are
# absent from result/ and listed with a reason in result/manifest.tsv and
# result/qc_status.csv; result/metrics.csv has one row per sample, with NA
# where a step did not run; both reports render.

set -euo pipefail

dryrun=()
if [[ ${1:-} == -n ]]; then dryrun=(-n); shift; fi
[[ $# -ge 2 ]] || { sed -n '6,13p' "$0" | sed 's/^# \{0,1\}//' >&2; exit 1; }
sheet=$1
outdir=$2
threads=${3:-32}
mkdir -p "$outdir"
outdir=$(cd "$outdir" && pwd)

awk 'BEGIN {FS=OFS="\t"}
	NR == 1 {
		for (i = 1; i <= NF; i++) col[$i] = i
		s = ("sample" in col) ? col["sample"] : col["sampleid"]
		c = ("cellranger" in col) ? col["cellranger"] : col["cellranger_dir"]
		if (!s || !c) {print "sample sheet needs sample/cellranger columns" > "/dev/stderr"; exit 1}
		print "sample", "cellranger"; next
	}
	NF {print $s, $c}' "$sheet" > "$outdir/samples.txt"
echo "$(($(wc -l < "$outdir/samples.txt") - 1)) samples -> $outdir/samples.txt"

# snRNA-seq settings: a stricter mito cut-off than the default, and the
# GEM-X 3' v4 multiplet rate (~0.4% per 1,000 nuclei recovered).
cat > "$outdir/config.yaml" <<'EOF'
seed: 42
ambient:
  method: soupx
  compare: [decontx]
filterbycount:
  mincount: 500
  minfeature: 300
  mito: 5
doublet:
  run: [doubletfinder, scdblfinder]
  decider: doubletfinder
  findpK: false
  pK: 0.01
  numthreads: 5
  rate: 0.08
  capacity: 20000
  nreaction: 1
EOF

cellqc -d "$outdir" -t "$threads" "${dryrun[@]}" -c "$outdir/config.yaml" -- "$outdir/samples.txt"
