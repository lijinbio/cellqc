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

