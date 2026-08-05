#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:
#
# Regenerate the figures in docs/. Run it after changing the workflow.
#
#   bash docs/make_figures.sh
#
# docs/workflow.png  -- hand-drawn pipeline overview, from docs/workflow.dot
# docs/tests/dag.png -- Snakemake job DAG, from the workflow itself, so it
#                       cannot describe a pipeline that no longer exists. It is
#                       built against empty stub Cell Ranger directories: --dag
#                       only needs the input paths to exist, not to contain
#                       anything.
#
# The DAG is drawn for ONE sample. Both report rules take every sample's outputs
# as input, so a two-sample DAG has ~40 edges converging on two nodes and reads
# as spaghetti; the per-sample rules simply repeat, which the caption says.
#
# Needs graphviz (dot) and, for the DAG, the cellqc environment on PATH.

set -euo pipefail

DPI=300
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo=$(dirname "$here")

command -v dot >/dev/null || { echo "graphviz 'dot' not found" >&2; exit 1; }

echo "==> docs/workflow.png"
dot -Tpng -Gdpi=$DPI "$here/workflow.dot" -o "$here/workflow.png"

if ! command -v snakemake >/dev/null; then
	echo "snakemake not found: skipping docs/tests/dag.png (activate the cellqc env)" >&2
	exit 0
fi

echo "==> docs/tests/dag.png"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

printf 'sample\tcellranger\tnreaction\n' > "$tmp/samples.txt"
for s in AMD1; do
	mkdir -p "$tmp/cellranger/$s/outs"
	touch "$tmp/cellranger/$s/outs/raw_feature_bc_matrix.h5" \
		"$tmp/cellranger/$s/outs/filtered_feature_bc_matrix.h5" \
		"$tmp/cellranger/$s/outs/possorted_genome_bam.bam" \
		"$tmp/cellranger/$s/outs/possorted_genome_bam.bam.bai"
	printf '%s\tcellranger/%s/outs\t1\n' "$s" "$s" >> "$tmp/samples.txt"
done

cat > "$tmp/config.yaml" <<'YAML'
seed: 42
ambient:
  method: soupx
  compare: [decontx]
doublet:
  run: [doubletfinder, scdblfinder]
  decider: doubletfinder
YAML

# config.smk prints the effective config to stdout before the graph, so keep
# only the digraph. concentrate=true merges the edge bundles converging on the
# two report rules, which is what makes the picture readable.
snakemake --snakefile "$repo/cellqc/Snakefile" --directory "$tmp" \
	--config samplefile="$tmp/samples.txt" outdir="$tmp/out" \
		configfile="$tmp/config.yaml" nowtimestr=docs \
	--configfile "$tmp/config.yaml" --dag 2>/dev/null \
	| sed -n '/^digraph/,$p' \
	| sed 's/graph\[bgcolor=white, margin=0\];/graph[bgcolor=white, margin=0.1, concentrate=true, nodesep=0.35, ranksep=0.55];/' \
	> "$tmp/dag.dot"

mkdir -p "$here/tests"
dot -Tpng -Gdpi=$DPI "$tmp/dag.dot" -o "$here/tests/dag.png"

echo "done."
