#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:
#
# Smoke test for an installed cellqc: builds the DAG for a stub cohort and
# checks the outputs the pipeline promises. No data, no cluster, seconds.
#
#   bash tests/dryrun.sh
#
# --dry-run only needs the input paths to exist, so empty files are enough.

set -uo pipefail

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
fail=0
say() { if [[ $1 -eq 0 ]]; then echo "  ok    $2"; else echo "  FAIL  $2"; fail=1; fi; }

for s in WITHBAM NOBAM; do
	mkdir -p "$tmp/cr/$s/outs"
	touch "$tmp/cr/$s/outs/raw_feature_bc_matrix.h5" "$tmp/cr/$s/outs/filtered_feature_bc_matrix.h5"
done
touch "$tmp/cr/WITHBAM/outs/possorted_genome_bam.bam" "$tmp/cr/WITHBAM/outs/possorted_genome_bam.bam.bai"

{
	printf 'sample\tcellranger\tnreaction\n'
	printf 'WITHBAM\tcr/WITHBAM/outs\t1\n'
	printf 'NOBAM\tcr/NOBAM/outs\t2\n'
} > "$tmp/samples.txt"

dry() { cellqc -d "$tmp/out" -t 2 -n "$@" -- "$tmp/samples.txt" > "$tmp/log" 2>&1; }

dry
say $? 'the workflow builds a DAG'

for f in result/WITHBAM.h5ad result/WITHBAM_obs.txt.gz result/WITHBAM_var.txt.gz \
	filterdoublet/WITHBAM.h5ad result/report.html result/report_slides.pdf result/metrics.csv; do
	grep -q "$f" "$tmp/log"
	say $? "produces $f"
done

[[ $(grep -c '^rule nuclear_fraction:' "$tmp/log") -eq 1 ]]
say $? 'nuclear fraction runs only for the sample with a BAM'

printf 'doublet:\n  skip: true\n' > "$tmp/removed.yaml"
dry -c "$tmp/removed.yaml"
[[ $? -ne 0 ]]
say $? 'a removed config key is rejected instead of ignored'

python -c 'from cellqc.qcutil import expected_doublet_rate as e; \
	assert e(13000, 2, 0.1, 13000)[0] == e(13000, 1, 0.1, 13000)[0] / 2'
say $? 'nreaction divides the expected doublet rate'

if [[ $fail -eq 0 ]]; then echo 'dryrun: PASS'; else echo 'dryrun: FAIL'; fi
exit $fail
