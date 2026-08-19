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

# No nreaction column at all: the value comes from doublet.nreaction.
{
	printf 'sample\tcellranger\n'
	printf 'WITHBAM\tcr/WITHBAM/outs\n'
} > "$tmp/samples_nonreaction.txt"

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

# The gene sets exist twice by necessity -- config.smk cannot import the package,
# and the effective definition has to reach config_<timestamp>.yaml. Compare the
# dumped config against qcutil so the two copies cannot drift silently.
python - "$tmp/out" <<'PY'
import glob, sys, yaml
from cellqc.qcutil import GENE_SETS
dumped = sorted(glob.glob(sys.argv[1] + '/config_*.yaml'))[-1]
cfg = yaml.safe_load(open(dumped))['geneset']
assert cfg == {k: dict(v) for k, v in GENE_SETS.items()}, (
	f'config.smk defaults and qcutil.GENE_SETS disagree:\n{cfg}\n{GENE_SETS}')
PY
say $? 'config.smk defaults and qcutil.GENE_SETS agree'

printf 'doublet:\n  skip: true\n' > "$tmp/removed.yaml"
dry -c "$tmp/removed.yaml"
[[ $? -ne 0 ]]
say $? 'a removed config key is rejected instead of ignored'

python -c 'from cellqc.qcutil import expected_doublet_rate as e; \
	assert e(13000, 2, 0.1, 13000)[0] == e(13000, 1, 0.1, 13000)[0] / 2'
say $? 'nreaction divides the expected doublet rate'

# Defining one set leaves the others at their defaults, so `mt` survives a
# config that only redefines `ribo` -- but emptying it is an error, because
# filterbycount.mito is a threshold on a metric that would no longer exist.
printf 'geneset:\n  ribo:\n    patterns: ["^RP[SL]"]\n' > "$tmp/oneset.yaml"
dry -c "$tmp/oneset.yaml"
say $? 'redefining one gene set keeps the other defaults'

printf 'geneset:\n  mt:\n' > "$tmp/nomt.yaml"
dry -c "$tmp/nomt.yaml"
[[ $? -ne 0 ]]
say $? 'emptying the mt gene set is rejected'

cellqc -d "$tmp/out_nr" -t 2 -n -- "$tmp/samples_nonreaction.txt" > "$tmp/log_nr" 2>&1
say $? 'a sample file with no nreaction column builds a DAG (doublet.nreaction supplies it)'

# The gene sets are what makes a non-human reference work, so assert the three
# naming conventions rather than only that the DAG builds: human/mouse match by
# prefix pattern, macaque falls back to bare mtDNA symbols, and the near-misses
# that share a prefix with a set stay out.
python - <<'PY'
from cellqc.qcutil import GENE_SETS, gene_set_mask
def hits(names, s):
	m, how = gene_set_mask(names, GENE_SETS[s])
	return sorted(n for n, k in zip(names, m) if k), how
human = ['MT-ND1', 'MT-CO1', 'RPSA', 'RPSA-1', 'RPS6', 'RPS6KA1', 'RPS19BP1',
	'HBB', 'HBA1', 'HBEGF', 'HBP1', 'HBS1L', 'ACTB']
mouse = ['mt-Nd1', 'Rpl13a', 'Hba-a1', 'Hbb-bs', 'Hbegf', 'Actb']
macaque = ['ND1', 'COX1', 'CYTB', 'PTGS1', 'MT2A', 'ACTB']
assert hits(human, 'mt') == (['MT-CO1', 'MT-ND1'], 'pattern')
assert hits(human, 'ribo') == (['RPS6', 'RPSA', 'RPSA-1'], 'pattern')
assert hits(human, 'hb') == (['HBA1', 'HBB'], 'pattern')
assert hits(mouse, 'mt') == (['mt-Nd1'], 'pattern')
assert hits(mouse, 'hb') == (['Hba-a1', 'Hbb-bs'], 'pattern')
assert hits(macaque, 'mt') == (['COX1', 'CYTB', 'ND1'], 'symbol')
# The fallback is a fallback: a reference with MT- genes never claims bare COX1.
assert hits(human + ['COX1'], 'mt') == (['MT-CO1', 'MT-ND1'], 'pattern')
PY
say $? 'gene sets match human, mouse and macaque and exclude the near-misses'

if [[ $fail -eq 0 ]]; then echo 'dryrun: PASS'; else echo 'dryrun: FAIL'; fi
exit $fail
