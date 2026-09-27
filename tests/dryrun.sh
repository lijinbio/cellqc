#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:
#
# Smoke test: build the DAG for two stub samples and check what it promises.
# No data needed -- a dry run only needs the input paths to exist.
#
#   bash tests/dryrun.sh

set -uo pipefail
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
fail=0
check() { if eval "$2" >/dev/null 2>&1; then echo "  ok    $1"; else echo "  FAIL  $1"; fail=1; fi; }

# Two stub samples, one with a BAM.
for s in WITHBAM NOBAM; do
	mkdir -p "$tmp/cr/$s"
	touch "$tmp/cr/$s/raw_feature_bc_matrix.h5" "$tmp/cr/$s/filtered_feature_bc_matrix.h5"
done
touch "$tmp/cr/WITHBAM/possorted_genome_bam.bam" "$tmp/cr/WITHBAM/possorted_genome_bam.bam.bai"
printf 'sample\tcellranger\nWITHBAM\tcr/WITHBAM\nNOBAM\tcr/NOBAM\n' > "$tmp/samples.txt"

dry() { cellqc -d "$tmp/out" -t 2 -n "$@" -- "$tmp/samples.txt" > "$tmp/log" 2>&1; }

check 'the DAG builds' 'dry'
for f in result/qc_status.csv result/manifest.tsv result/metrics.csv result/report.html result/report_slides.pdf; do
	check "produces $f" "grep -q $f $tmp/log"
done
check 'every step writes a status file' \
	'for s in ambient filterbycount doubletfinder scdblfinder filterdoublet postproc; do grep -q $s/WITHBAM_status.tsv $tmp/log || exit 1; done'
check 'nuclear fraction runs only where there is a BAM' '[[ $(grep -c "^rule nuclear_fraction:" $tmp/log) -eq 1 ]]'

# config.smk and qcutil.GENE_SETS hold the same defaults; they must not drift.
check 'config.smk and qcutil.GENE_SETS agree' "python -c '
import glob, yaml; from cellqc.qcutil import GENE_SETS
cfg = yaml.safe_load(open(sorted(glob.glob(\"$tmp/out/config_*.yaml\"))[-1]))
assert cfg[\"geneset\"] == {k: dict(v) for k, v in GENE_SETS.items()}'"

printf 'doublet:\n  skip: true\n' > "$tmp/bad.yaml"
check 'a removed config key is rejected' '! dry -c $tmp/bad.yaml'

[[ $fail -eq 0 ]] && echo 'dryrun: PASS' || echo 'dryrun: FAIL'
exit $fail
