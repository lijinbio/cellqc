#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:
#
# Data-free check of the config layer and the DAG shape. Runs in seconds, needs
# no Cell Ranger output, no HPC and no test data: `--dag`/`--dry-run` only needs
# the input paths to exist, so stub files are enough.
#
#   bash tests/dryrun.sh
#
# What it covers is exactly what an end-to-end run is too slow to iterate on:
# which rules the DAG contains for a given config, and whether a stale or
# invalid config is rejected with a message that says what to do. The numerical
# behaviour of the steps is not tested here -- that is the reference run in
# tests/mwe/.

set -uo pipefail

pass=0 fail=0
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

command -v cellqc >/dev/null || { echo "cellqc not on PATH: activate the cellqc environment" >&2; exit 1; }

# Two stub Cell Ranger directories; only WITHBAM carries an indexed BAM, so the
# same run exercises both branches of the nuclear-fraction detection.
mk_cellranger() { # <dir> <bam:yes|no>
	mkdir -p "$1"
	touch "$1/raw_feature_bc_matrix.h5" "$1/filtered_feature_bc_matrix.h5"
	if [[ $2 == yes ]]; then
		touch "$1/possorted_genome_bam.bam" "$1/possorted_genome_bam.bam.bai"
	fi
}
mk_cellranger "$tmp/cr/WITHBAM/outs" yes
mk_cellranger "$tmp/cr/NOBAM/outs" no

printf 'sample\tcellranger\tnreaction\n' > "$tmp/samples.txt"
printf 'WITHBAM\tcr/WITHBAM/outs\t1\n' >> "$tmp/samples.txt"
printf 'NOBAM\tcr/NOBAM/outs\t2\n' >> "$tmp/samples.txt"

n=0
dryrun() { # <config file or ''>  -> writes combined output to $out
	n=$((n + 1))
	out="$tmp/run$n.log"
	if [[ -n ${1:-} ]]; then
		cellqc -d "$tmp/out$n" -t 2 -n -c "$1" -- "$tmp/samples.txt" > "$out" 2>&1
	else
		cellqc -d "$tmp/out$n" -t 2 -n -- "$tmp/samples.txt" > "$out" 2>&1
	fi
	return $?
}

ok()   { pass=$((pass + 1)); printf '  ok    %s\n' "$1"; }
bad()  { fail=$((fail + 1)); printf '  FAIL  %s\n' "$1"; [[ -n ${2:-} ]] && sed -n '1,40p' "$2"; }
check() { # <description> <condition-exit-code> [logfile]
	if [[ $2 -eq 0 ]]; then ok "$1"; else bad "$1" "${3:-}"; fi
}

echo '== defaults: every stage present, nuclear fraction only where there is a BAM =='
dryrun ''
rc=$?
check 'dry run succeeds' $rc "$out"
for rule in ambient barcoderank filterbycount doubletfinder scdblfinder filterdoublet postproc qcreport slidereport; do
	grep -qE "^rule $rule:" "$out"; check "rule $rule in the DAG" $? "$out"
done
# One nuclear_fraction job, for WITHBAM only.
[[ $(grep -cE '^rule nuclear_fraction:' "$out") -eq 1 ]]
check 'nuclear_fraction runs once (BAM present for one sample)' $? "$out"
grep -q 'no indexed possorted_genome_bam.bam' "$out"
check 'the sample without a BAM is reported, not silently skipped' $? "$out"

echo '== doublet.run selects the callers, with no skip flag =='
cat > "$tmp/one_caller.yaml" <<'YAML'
doublet:
  run: [scdblfinder]
  decider: scdblfinder
YAML
dryrun "$tmp/one_caller.yaml"
check 'dry run succeeds with a single caller' $? "$out"
grep -qE '^rule scdblfinder:' "$out"; check 'scdblfinder runs' $? "$out"
! grep -qE '^rule doubletfinder:' "$out"; check 'doubletfinder does not run' $? "$out"

echo '== stale configs =='
cat > "$tmp/legacy.yaml" <<'YAML'
dropkick:
  skip: true
scpred:
  reference: /nonexistent/model.rds
doubletfinder:
  pK: 0.02
  findpK: false
YAML
dryrun "$tmp/legacy.yaml"
check 'a v0.1.0 config still runs' $? "$out"
grep -q "section 'dropkick' is obsolete" "$out"; check 'dropkick warns' $? "$out"
grep -q "section 'scpred' is obsolete" "$out"; check 'scpred warns' $? "$out"
grep -q "'doubletfinder' was renamed to 'doublet'" "$out"; check 'doubletfinder section migrates' $? "$out"
grep -q '"pK": 0.02' "$out"; check 'the migrated pK is kept' $? "$out"

printf 'doublet:\n  skip: false\n' > "$tmp/skipfalse.yaml"
dryrun "$tmp/skipfalse.yaml"
check "'doublet.skip: false' still runs" $? "$out"
grep -q "'doublet.skip' is obsolete" "$out"; check "'doublet.skip: false' warns" $? "$out"

echo '== configs that must be rejected =='
printf 'doublet:\n  skip: true\n' > "$tmp/skiptrue.yaml"
dryrun "$tmp/skiptrue.yaml"
[[ $? -ne 0 ]]; check "'doublet.skip: true' is an error, not a silent behaviour change" $?
grep -q "doublet detection always runs" "$out"; check 'the error says what to do instead' $? "$out"

cat > "$tmp/baddecider.yaml" <<'YAML'
doublet:
  run: [scdblfinder]
  decider: doubletfinder
YAML
dryrun "$tmp/baddecider.yaml"
[[ $? -ne 0 ]]; check 'a decider that is not in run is rejected' $?

printf 'filterbycount:\n  mito: 200\n' > "$tmp/badmito.yaml"
dryrun "$tmp/badmito.yaml"
[[ $? -ne 0 ]]; check 'an out-of-range threshold is rejected by the schema' $?

echo
echo "passed $pass, failed $fail"
[[ $fail -eq 0 ]]
