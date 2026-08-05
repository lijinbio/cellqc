# v0.3.0 - Aug 5, 2026

## Breaking changes

- Output layout. The final matrix is now `result/{sample}.h5ad` — postproc's output, what a user actually
  takes away — and the pre-integration matrix moved from `result/` to `filterdoublet/{sample}.h5ad`, where
  every other stage's output already lives. The `postproc/` directory is gone. Downstream code that read
  `postproc/{sample}.h5ad` should read `result/{sample}.h5ad`; code that read the old `result/{sample}.h5ad`
  wants `filterdoublet/{sample}.h5ad`, which is now `temp()` — see below. This resolves the open question
  in `docs/v0.2.0_plan.md` §9b.3.
- Removed `doublet.skip`. Doublet detection always runs; the callers are `doublet.run` and a caller you do
  not want is left out of it, which is the same convention `nuclear_fraction` already used (no skip flag).
  `doublet.skip: false` warns and is dropped. `doublet.skip: true` is an **error**, not a warning: there is
  no configuration that keeps every called doublet, so continuing would quietly remove cells from a run
  that asked for none.

## Changed

- Packaging moved from `setup.py` to a PEP 621 `pyproject.toml` (setuptools backend). Metadata is
  declarative, the version is still single-sourced from `cellqc/__init__.py`, the license is an SPDX
  expression, and the workflow files ship via `[tool.setuptools.package-data]`. `MANIFEST.in` now only adds
  what an sdist needs beyond the package, and no longer sweeps `CLAUDE.md` into the distribution.
- The nuclear fraction stays on the pysam implementation. Switching to DropletQC was considered and
  rejected: it is unmaintained (still `0.0.0.9000`, never released), it is GitHub-only, and its dependency
  closure adds ~15 R packages (`GenomicFeatures`, `rtracklayer`, `ggpubr`, …) to an environment that does
  not otherwise need them. `tests/validate_nuclear_fraction.py` remains as the evidence that the two agree
  (r = 1.000000, max |Δ| = 5.7e-4 over 13,559 barcodes).
- Figures: PNG output raised from 200 to 300 dpi and rasterized layers inside the vector PDFs from 500 to
  600 dpi. Both reports get the higher-resolution PNGs, so `result/report.html` grows accordingly. The R
  steps (`ggsave`, `png`) carry the same numbers — change them together with `cellqc/qcutil.py`.
- Slide deck: the "Cells retained at each stage" table ran off the right-hand edge of the slide and the
  Cell Ranger metrics were set in 5pt type in the middle of an empty frame. Tables now get real column
  names, are folded two-up where they are long, and are wrapped in a fit-to-width box that shrinks a table
  only when it would overflow, so no data-dependent table can silently run off a slide again.
- `docs/workflow.png` redrawn (it still showed dropkick and scPred) and now generated from
  `docs/workflow.dot` by `bash docs/make_figures.sh`. The Snakemake job-DAG image is gone: it duplicated
  what the workflow diagram already shows, at rule-name granularity nobody reads, and went stale on every
  rule rename. `docs/tests/` carries the current example report, slide deck and metrics.csv.

## Added

- `result/metrics.csv`: every scalar the run produced, one row per sample — Cell Ranger metrics, knee and
  inflection, ambient contamination per method, per-criterion filter counts, each doublet caller's count
  and their concordance, nuclear-fraction quartiles, retained fraction. Assembled from the same
  `reportdata.collect()` the reports use, so it cannot disagree with them, and it means nothing downstream
  has to scrape a number out of the HTML. 76 columns on the reference sample.
- `.obs` and `.var` are written beside `result/{sample}.h5ad` as gzipped TSVs (`{sample}_obs.txt.gz`,
  `{sample}_var.txt.gz`), indexed by `barcode` and `gene`, so the per-cell QC metrics and the feature table
  can be read and joined without anndata.
- `filterdoublet/{sample}.h5ad` is `temp()`: Snakemake deletes it once `result/{sample}.h5ad` is written.
  It held the same cells and the same counts as the final matrix, differing only in the barcode prefix,
  the uniquified var names and the nuclear-fraction columns, so keeping it wrote every count matrix to
  disk twice — 72 MB per sample as written on the reference cohort (~35 MB of actual disk there, since
  that filesystem compresses), against ~207 MB for the whole run. `snakemake --notemp` keeps it.
- README documents *why* the expected doublet rate is linear in cell yield — droplet occupancy is Poisson,
  so the multiplet fraction among occupied droplets is `1 - λ/(e^λ - 1) ≈ λ/2` — with the references
  behind the rule of thumb (Bloom 2018 PeerJ 6:e5578; 10x Chromium user guides, ≈0.8% per 1,000 cells
  recovered; scDblFinder's ≈1% per 1,000) and the two known upward biases.
- `tests/dryrun.sh`: a data-free smoke test — the DAG builds, the promised outputs are produced, the
  nuclear fraction runs only where there is a BAM, an obsolete config key is rejected, and `nreaction`
  scales the expected doublet rate. Seconds, no HPC, no test data.

## Removed

- `tests/` is now two scripts and a README, with no lab-specific paths, accounts or helpers, so it can ship
  publicly and is purely about testing the package. Gone: `tests/main.sh` (depended on `trapdebug`,
  `mrrdir.sh`, `slurmtaco.sh` and called `cellqc` without a sample file), `tests/mwe/slurm_cellqc.sh` (a
  site-specific sbatch submission; the reference run is documented in `tests/README.md` instead) and
  `tests/nreaction/` (its one assertion moved into `dryrun.sh`).
  `tests/mwe/validate_nuclear_fraction.py` moved to `tests/validate_nuclear_fraction.py`.

# v0.2.0 - Aug 5, 2026

## Breaking changes

- Removed dropkick (empty-droplet calling) and scPred (cell-type annotation). Cell calling is now Cell
  Ranger EmptyDrops alone; annotate downstream of cellqc. Old configs carrying `dropkick:`/`scpred:`
  sections warn and continue rather than failing.
- Removed SeuratDisk. `.h5ad` is the only interchange format; `.h5seurat` intermediates and
  `result/*.h5seurat` are gone. This lifts the Seurat v4 pin -- v0.2.0 targets Seurat 5.
- The `doubletfinder:` config section is renamed `doublet:` (old configs are migrated with a warning),
  because the step now supports more than one caller.
- `cellqc/cellqc.py` renamed to `cellqc/cli.py` (entry point `cellqc.cli:main`). A module sharing the
  package name shadowed the package whenever Snakemake put the workflow directory on `sys.path`.

## Reproducibility (important)

- A random seed is now set for every stochastic step (`seed`, default 42). v0.1.0 set no seed anywhere.
  This was worse than a cosmetic issue: `SoupX::adjustCounts(roundToInt=TRUE)` uses randomised rounding,
  so **every v0.1.0 run produced a slightly different integer count matrix**, which propagated into the
  mitochondrial percentage and moved cells across the filter threshold. Verified: with the seed fixed,
  two runs are identical; changing the seed shifts the corrected total by ~1,300 of 91.3M counts.
- Consequence: v0.2.0 results are not bit-for-bit comparable with v0.1.0. On the reference sample the
  final cell count differs by 1 (10,223 vs 10,222) for exactly this reason.

## Added

- PDF slide report (`result/report_slides.pdf`), beamer via tectonic, generated from a data-driven Jinja2
  template so new samples need no template edit. Includes the Cell Ranger metrics table and a barcode
  rank plot.
- Barcode rank (knee) plot with knee/inflection computed the DropletUtils way. Diagnostic only.
- Nuclear fraction vs log10(UMI) scatter plot.
- DecontX as an alternative ambient method (`ambient.method`), plus `ambient.compare` to report other
  methods' contamination estimates without applying them.
- scDblFinder as a second doublet caller. `doublet.run` selects which callers execute and
  `doublet.decider` selects the single caller that removes cells; every caller's score lands in `.obs`.
  Caller concordance (2x2 table and Cohen's kappa) is reported.
- Per-criterion filtering counts: how many cells fail each of mincount/minfeature/mito individually,
  only, and in combination. v0.1.0 reported only the before/after totals.
- Ambient correction impact (counts removed) is now reported; v0.1.0 applied it silently.
- Both reports carry an explicit limitations section.
- All figures are emitted as vector PDF with editable text (dense layers rasterized at 500 dpi) plus PNG
  for the HTML report.
- `envs/cellqc.yaml`: a single conda environment for the whole pipeline.

## Changed

- The nuclear fraction is reimplemented in pysam, removing the GitHub-only DropletQC dependency.
  Validated against DropletQC on the reference sample: Pearson r = 1.000000, median |delta| = 0.000000,
  max |delta| = 0.00057 across 13,559 barcodes.
- The nuclear-fraction step is enabled automatically per sample by the presence of an indexed
  `possorted_genome_bam.bam`. There is no skip flag, and cohorts with mixed BAM availability work.
- `filterbycount` moved from Seurat to scanpy. The criteria and the `^MT-|^mt-` pattern are unchanged, so
  thresholds carry over.
- The expected-doublet-rate constants (0.1 and 13000) are exposed as `doublet.rate` and
  `doublet.capacity` instead of being hard-coded. Defaults reproduce v0.1.0.
- The config schema now validates types and ranges for every parameter instead of only requiring
  `samples`.

## Fixed

- `doubletFinder(reuse.pANN=FALSE)` now passes `NULL`. Upstream DoubletFinder v2.0.6 changed this check
  from `if (reuse.pANN)` to `if (!is.null(reuse.pANN))`, so `FALSE` took the reuse branch and failed with
  "cannot xtfrm data frames".
- The `lijinbio/DoubletFinder` fork is no longer needed: upstream v2.0.6 handles Seurat 5 layers.
- Ensembl gene IDs are preserved through the ambient step. v0.1.0 went through `Seurat::Read10X`, which
  keeps only made-unique gene symbols.

# v0.0.9 - Apr 5, 2025

- Calculate a nuclear fraction score to quantify the proportion of reads from intronic regions.

# v0.0.8 - Dec 5, 2024

- Add the sample ID to the cell barcode prefix in the postprocessed .h5ad file.

# v0.0.7 - Jan 29, 2024

- Bug fix: Updated DoubletFinder with modified function names.
- Added "-D|--define" to support defining an individual sample without a SAMPLEFILE.

# v0.0.6 - Mar 8, 2023

- Added support for "skip=True" in DoubletFinder, useful for multiplex libraries.

# v0.0.5 - Feb 16, 2023

- Separated sample file (e.g., samples.txt) from config.yaml.

# v0.0.4 - Dec 16, 2022

- Updated installation instructions using conda.

# v0.0.3 - Oct 31, 2022

## Initial Implementation

- Implemented conditional execution for Dropkick and scPred.
- Implemented qc_report.html for a QC summary.

