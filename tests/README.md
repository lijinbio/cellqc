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
