# tests

There is no unit-test framework. Testing is three things: a fast structural check that needs no data, an
end-to-end run on a reference sample, and a hard numerical gate on the one component that replaced an
external implementation.

## Structural check (seconds, no data)

```bash
bash tests/dryrun.sh
```

Builds stub Cell Ranger directories — `--dry-run` only needs the input paths to exist — and asserts what
an end-to-end run is too slow to iterate on:

- every stage appears in the DAG, and `nuclear_fraction` appears only for the sample with an indexed BAM;
- `doublet.run` selects the callers (there is no skip flag);
- v0.1.0 configs still run: `dropkick`/`scpred` warn and are dropped, `doubletfinder:` migrates to
  `doublet:` keeping its values;
- configs that must be rejected are rejected: `doublet.skip: true` (removed — honouring it silently is
  impossible now), a `decider` that is not in `run`, and an out-of-range threshold.

Prints a pass/fail count and exits non-zero if anything failed.

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

Expected on this sample: 13,559 → 11,234 (count filter) → 10,223 cells; DoubletFinder 1,011 (9.00%) vs
scDblFinder 1,153 (10.26%), Cohen's κ = 0.759. The same seed must give an identical matrix.

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

- `nreaction/main.sh` — `nreaction > 1`, which divides the expected doublet rate per reaction. Checks the
  arithmetic without any data; give it `CELLRANGER=/path/to/outs` to run the pipeline as well.

## Example outputs

`docs/tests/` holds the HTML report and the slide deck from the reference run, plus the job DAG.
Regenerate the DAG and the workflow diagram with `bash docs/make_figures.sh`.
