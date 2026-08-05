# tests

There is no unit-test framework. Validation is an end-to-end run on a reference sample plus a hard
numerical gate on the one component that replaced an external implementation.

## End-to-end run

`mwe/slurm_cellqc.sh` runs the full pipeline on GSE188280_GSM5676874_0715_Macula_Retina
(Cell Ranger 10.0.0, 13,559 called cells) and keeps the exact submission for the record.

```bash
WORK=/path/on/shared/storage        # NOT /tmp: it is node-local and invisible to compute nodes
mkdir -p "$WORK"
sbatch --chdir="$WORK" tests/mwe/slurm_cellqc.sh
```

Free partition, `ruic20_lab` account, `--requeue` so preemption resubmits. Snakemake resumes from
completed outputs, so a requeue only costs the interrupted rule.

## Nuclear-fraction acceptance gate

v0.2.0 replaced `DropletQC::nuclear_fraction_tags()` with a pysam implementation to drop a GitHub-only
dependency. That is only acceptable if it reproduces the reference:

```bash
python tests/mwe/validate_nuclear_fraction.py \
  "$WORK"/out/nuclear_fraction/<sample>.txt.gz \
  <v0.1.0 DropletQC output>.txt.gz \
  "$WORK"
```

Gate: identical barcode sets, Pearson and Spearman > 0.999, median |Δ| < 0.001, max |Δ| < 0.01.
**If it fails, revert to DropletQC — do not loosen the thresholds.**

Last result: r = 1.000000, median |Δ| = 0.000000, max |Δ| = 0.000571 over 13,559 barcodes. Passed.

## Scenarios

- `nreaction/` — `nreaction > 1`, which lowers the expected doublet rate per reaction.

The `main.sh` scripts depend on lab HPC helpers (`trapdebug`, `mrrdir.sh`, `slurmtaco.sh`) that are not
in this repository.
