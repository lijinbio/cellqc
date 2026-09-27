# tests

Three scripts, no framework, no test data in the repository.

## `dryrun.sh` — smoke test

```bash
bash tests/dryrun.sh
```

Builds two stub Cell Ranger directories, one with a BAM, runs `cellqc -n` over them, and checks that:
- the DAG builds;
- the cohort outputs are promised (status tables, metrics, both reports);
- every per-sample step writes a status file;
- the nuclear fraction runs only for the sample with a BAM;
- the gene-set defaults in `config.smk` and `qcutil.GENE_SETS` agree;
- a removed config key is rejected.

Seconds, no data. Run it after changing `rules/`, the schema, `Snakefile` or `qcutil.GENE_SETS`.

## `main.sh` — failure tolerance on a real cohort

```bash
bash tests/main.sh -n samples.tsv out   # dry run
bash tests/main.sh samples.tsv out 32   # run; the same command resumes
```

Runs a whole cohort with snRNA-seq settings, which `main.sh` writes to `out/config.yaml`. The sample sheet
needs `sample` and `cellranger` columns; `sampleid` and `cellranger_dir` are accepted too. Use a cohort
that includes low-quality libraries: ones where SoupX finds no marker genes, or where filtering leaves too
few cells for doublet detection. Before v0.3.4 either failure stopped the run for every sample. Pass
criteria:
- the run exits 0;
- the failing samples are absent from `result/` and listed with a reason in `result/manifest.tsv` and
  `result/qc_status.csv`;
- `result/metrics.csv` has one row per sample, with `NA` for what was not computed;
- both reports render, with a *Sample status* section and a placeholder in place of each missing figure.

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
