# tests

Two scripts, no framework, no test data in the repository.

## `dryrun.sh` — smoke test

```bash
bash tests/dryrun.sh
```

Builds stub Cell Ranger directories, runs `cellqc -n` over them, and checks that the workflow still
produces what it promises: the final matrix and its `.obs`/`.var` dumps, the pre-integration matrix, both
reports, the nuclear fraction only for the sample that has a BAM, a rejected obsolete config key, and the
`nreaction` scaling of the expected doublet rate. Seconds, no cluster, no data — `--dry-run` only needs
the input paths to exist. Prints `PASS`/`FAIL` and exits non-zero on failure.

Run it after changing `rules/config.smk`, the schema, `Snakefile`, or any rule's outputs.

## `validate_nuclear_fraction.py` — acceptance gate

cellqc computes the nuclear fraction with pysam rather than depending on
[DropletQC](https://github.com/powellgenomicslab/DropletQC). That is only acceptable if it reproduces the
reference implementation, so this compares the two on the same sample:

```bash
python tests/validate_nuclear_fraction.py <cellqc>.txt.gz <dropletqc>.txt.gz [outdir]
```

Gate: identical barcode sets, Pearson and Spearman > 0.999, median |Δ| < 0.001, max |Δ| < 0.01. With an
`outdir` it also writes an agreement and Bland–Altman plot. **If it fails, the fix is to revert to
DropletQC, not to loosen the thresholds.**

Last run on GSE188280_GSM5676874_0715_Macula_Retina (Cell Ranger 10.0.0, 13,559 barcodes): r = 1.000000,
median |Δ| = 0.000000, max |Δ| = 0.000571. Passed.

## Reference run

The numbers quoted in the README and in `CHANGELOG.md` come from an end-to-end run on
[GSE188280](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188280) (sample GSM5676874, macula
retina), which anyone can reproduce:

```bash
printf 'sample\tcellranger\tnreaction\n' > samples.txt
printf 'GSE188280_Macula\t/path/to/cellranger/outs\t1\n' >> samples.txt
cellqc -d out -t 16 -- samples.txt
```

Expected: 13,559 → 11,234 cells (count filter) → 10,223 (doublets); DoubletFinder 1,011 (9.00%) vs
scDblFinder 1,153 (10.26%), Cohen's κ = 0.759. With the default seed the run is reproducible — a repeat
gives an identical matrix. Example outputs are in `docs/tests/`.

### Resuming an interrupted run

On a preemptible queue this run *will* be interrupted — the BAM pass alone is minutes of the wall clock.
Snakemake resumes from completed outputs, so a restart only costs the rule that was in flight, but two
things get in the way and neither is obvious from the error:

1. **A killed run leaves a stale lock.** Every restart then dies with `LockException: Directory cannot be
   locked` before doing any work. Clear it with `--unlock` first.
2. **The interrupted rule left a half-written output.** Snakemake refuses to trust it and asks for
   `--rerun-incomplete`.

The `cellqc` CLI does not pass Snakemake flags through, so a resume calls Snakemake directly with the
arguments the CLI builds:

```bash
pkg=$(python -c "import cellqc, pathlib; print(pathlib.Path(cellqc.__file__).parent)")
common=(--snakefile "$pkg/Snakefile" --directory out
        --config samplefile=samples.txt outdir=out configfile=config.yaml nowtimestr=resume
        --configfile config.yaml)

snakemake "${common[@]}" --unlock                       # only needed after a kill
snakemake "${common[@]}" --cores 16 --jobs 16 --rerun-incomplete
```

Do not delete the output directory to "start clean" — that throws away every completed stage and makes the
next preemption cost the whole run again.

Note also that `#SBATCH --requeue` is not a substitute: on the cluster this was developed on, preempted
jobs ended in state `PREEMPTED` and were never resubmitted, so the resume has to be driven by hand or by a
resubmit loop.
