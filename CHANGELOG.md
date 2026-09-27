# v0.3.5 - Sep 26, 2026

- No `postproc/` directory is left after a run. Its status file, `postproc/{sample}_status.tsv`, was the
  one `postproc` output not marked `temp`, so the directory survived holding only status files after the
  matrices were published to `result/`. The status reaches `result/qc_status.csv` through the `qcstatus`
  checkpoint, its only reader, and is now `temp` too. Snakemake removes the empty directory at the end of
  the run. Status files left in an existing outdir by v0.3.4 are not removed; delete them by hand.

# v0.3.4 - Sep 26, 2026

Per-sample failure tolerance. One sample failing one step no longer stops the run for every sample. A
sample for which every step succeeds is processed exactly as in v0.3.3, with the same matrix and the same
numbers. What changes is how a failure is handled, and a few output locations and columns, listed at the
end.

- **Every per-sample step records its outcome** in `{stage}/{sample}_status.tsv` with columns
  `sample, step, status, message`. `status` is `ok`, `fallback`, `failed` or `skipped`. One mechanism
  covers all steps, `cellqc/qcstatus.py` for the Python ones and its twin `scripts/guard.R` for the R ones,
  rather than a `try` in each script. The step body is guarded, but imports and `library()` calls are not,
  so a broken environment still stops the run. A failed step logs its traceback and its reason. A step whose
  input step is unusable is `skipped` and carries that step's reason forward, so the last step in a chain
  still names the root cause. Every declared output is left behind, as a 0-byte placeholder if the step
  did not write it, which keeps the DAG shape fixed. A matrix a failed step may have half-written is
  truncated rather than kept.
- **Fallbacks, where one exists.** If the applied ambient method fails, for example SoupX's `autoEstCont`
  with "No plausible marker genes found" on a low-complexity channel, the sample continues on the
  uncorrected Cell Ranger counts. This is recorded as `ambient: fallback`. The contamination table keeps a
  row for the failed method with `NA` estimates and adds a `none` row for what was applied. The ambient
  figure states the failure and its reason. A comparison-only method that fails is recorded as
  `ambient.<method>: failed` and changes nothing else.
- **Diagnostic steps never exclude a sample.** This covers `barcoderank` and `nuclear_fraction`: without
  the latter the final matrix is written without the nuclear fraction, as for a sample with no BAM. It also
  covers a non-deciding doublet caller: its `.obs` columns, summary row and concordance are `NA`, and the
  decider alone decides. A failed **decider** excludes the sample rather than handing the decision to the
  other caller, which would make one sample's doublet removal differ in method from the cohort's.
  Failures in `filterbycount` (including "all cells removed"), `filterdoublet` or `postproc` also exclude
  the sample.
- **`result/` holds only the samples that made it through.** A new `qcstatus` checkpoint collects every
  status into `result/qc_status.csv`, with every sample × step plus a final `result` row per sample, and
  into `result/manifest.tsv`, with columns `included`, `ncell`, `reason` and `fallback`. The final matrix is
  published to `result/` by a new `publish` rule only for included samples. It is a hard link to
  `postproc/`'s `temp()` output, so it is still written once, and no placeholder ever reaches `result/`.
- **Every sample in every table.** `result/metrics.csv` and every report table keep a row for each
  sample, with `NA` rather than an empty cell for what was not computed. `metrics.csv` gains
  `status_<step>`, `qc_included` and `qc_excluded_reason` as its last columns. A caller's doublet count now
  comes from the caller itself when `filterdoublet` did not run. The cell-retention cascade gains an
  `in_result` column.
- **Reports render from whatever is available.** Both reports open with a *Sample status* section: how
  many samples reached `result/`, the excluded ones with reasons, and every step that did not complete
  normally. The HTML report links to `qc_status.csv`. A figure a failed step did not produce keeps its slot
  with the reason, for example "Ambient RNA not available: ambient failed: ...", instead of vanishing. A
  missing, empty or unreadable stats file never raises. Messages are folded to ASCII for the slide deck, so
  an R error with a curly quote cannot break the LaTeX build.
- **The slide deck fits any cohort size.** The *Cells retained at each stage* table used to run off the
  bottom of its slide for about 15 or more samples, cutting off the table and the note under it. It is now
  split across frames of 12 rows. Every table in the deck is also scaled down when it is too tall, as it
  already was when too wide.
- **Exit status.** The run succeeds when at least one sample reaches `result/`. The excluded samples are
  logged at the end with their reasons. It fails when no sample does, and on any error that is not about
  one sample: a bad configuration, a missing Cell Ranger directory or matrix (now checked before anything
  runs), a missing package, or a deck that does not compile. The CLI passes `--keep-going` to Snakemake, so
  a job killed by the system does not stop the other samples' jobs.
- Output locations and columns. `result/{sample}_doublet_summary.txt` and `_doublet_concordance.txt` are
  written to `filterdoublet/` and published to `result/` with the matrix, so an excluded sample leaves no
  stub in `result/`. The final matrix's working copy is `postproc/{sample}.h5ad` (`temp`). A dry run
  lists `postproc/` rather than `result/` matrices, because which samples are published is decided at run
  time. `min_version` is 8.0, the floor `envs/cellqc.yaml` already declared.

# v0.3.3 - Aug 18, 2026

Reference-agnostic QC metrics. Additive for human and mouse: those runs are unchanged cell for cell,
`pct_counts_mt` being computed over the same genes as in v0.3.2 (two names change, listed at the end).
Macaque — and any other reference that names its mitochondrial genes without a contig prefix — now gets a
real mitochondrial percentage instead of a column of zeros.

- **QC gene sets are configuration** (`geneset` in the YAML), not a hard-coded pattern. A set is matched
  by case-insensitive `patterns`, with exact `symbols` as a fallback used only when no pattern matched
  anything, minus anything in `exclude`. Each set becomes a `pct_counts_<set>` column in `.obs`, a violin
  panel, and `n_<set>_genes` / `<set>_matched_by` / `median_pct_counts_<set>` in the filter statistics and
  `result/metrics.csv`. Adding a set is a config key, not a code change.
- **Macaque support.** `MT-`/`mt-` matches nothing in Ensembl Mmul_10, which names the 13 protein-coding
  mtDNA genes bare: `ND1`, `ND2`, `ND3`, `ND4`, `ND4L`, `ND5`, `ND6`, `COX1`, `COX2`, `COX3`, `ATP6`,
  `ATP8`, `CYTB`. Those are the default `geneset.mt.symbols`, tried only when the prefix pattern found
  nothing — a bare `COX1` is also a legacy alias of the nuclear gene *PTGS1*, so it is claimed only in a
  reference with no prefixed mitochondrial genes at all. Verified on the macaque retina reference: 13/13
  matched (`matched_by: symbol`), while GRCh38 and GRCm39 still match by pattern and never reach the
  fallback. Before this release such a run filtered on a `pct_counts_mt` that was 0 for every cell.
- **Ribosomal and hemoglobin percentages** are computed, plotted and recorded for every sample:
  `pct_counts_ribo` (`^RP[SL]\d`, `^RPLP\d`, `^RPSA$`, less `^RPS6K` and `^RPS19BP` — kinases and a
  binding protein, not ribosomal proteins) and `pct_counts_hb` (the globin cluster, written as full
  matches so `HBEGF`, `HBP1` and `HBS1L` are not swallowed by a bare `^HB`). **Neither is filtered on.** A
  defensible cut-off for either is tissue-dependent, and cellqc does not exclude cells on a number nobody
  has looked at; their violin panels carry no threshold line and are labelled *(not filtered)*.
- The QC violin figure grows a panel per gene set (five by default, was three) and widens with it instead
  of squeezing the existing panels. Both reports say in the figure caption that the panels without a
  dashed line remove no cells.
- Gene-set matching sees past `var_names_make_unique()`: a reference carrying `RPSA` twice becomes
  `RPSA` + `RPSA-1`, and an anchored pattern would otherwise count only the first copy.
- **`doublet.nreaction`** sets the per-sample default for every sample whose row in the sample file does
  not give its own; the column still wins where present. A cohort on one chemistry states it once — 10x
  GEM-X 3' v4 halves the multiplet rate of 3' v3.1 at the same yield, which `nreaction: 2` expresses
  against the unchanged `rate`/`capacity` line.
- **Packaging metadata is `pyproject.toml`'s alone.** The version moved out of `cellqc/__init__.py` to a
  static `[project] version`, and `__init__.py` reads it back with `importlib.metadata` so
  `cellqc.__version__`, `cellqc --version` and the reports are unchanged. `__author__`/`__email__` are
  gone — nothing read them, no PEP defines them, and `pyproject.toml`'s `authors` already carried the
  same values. Note for developers: a version bump now needs a `pip install -e .` before the reports show
  it, because `importlib.metadata` reads installed metadata rather than the working tree.
- `filter_n_mito_genes` in `result/metrics.csv` is now `filter_n_mt_genes`, one instance of the generic
  `n_<set>_genes` naming. `qcutil.MITO_PREFIXES` is replaced by `qcutil.GENE_SETS` and
  `qcutil.gene_set_mask()`; `qcutil.mito_percent()` returns `(pct, n_genes, matched_by)` rather than
  `(pct, n_genes)`.

# v0.3.2 - Aug 7, 2026

Additive release. Every cell that came out of v0.3.1 still comes out, with the same counts: the new `.obs`
columns are recorded, not filtered on, and the rest is what the reports show and say.

- Both reports call the tool **CellQC** — the HTML `<title>`, both report headings, the version lines and
  the prose that names it. The package, the command and the import stay lowercase `cellqc`; this is the
  product name as it appears to a reader, nothing else.
- The nuclear-fraction scatter (`nuclear_fraction/{sample}_nf_umi.{pdf,png}`) colours each cell by its
  mitochondrial percentage, which separates the two readings of the same corner of the plot: low depth
  with a high nuclear fraction is a damaged cell when mitochondrial content is high and a free nucleus or
  empty drop when it is not. The percentage comes from the same `filtered_feature_bc_matrix.h5` the plot
  already reads for total UMI, so it is the pre-ambient-correction number (`raw_pct_counts_mt` in the
  final `.obs`), and the colour scale is capped at the 99th percentile with an `extend='max'` arrow rather
  than letting a few near-100% cells flatten it. A reference with no gene matching `MT-`/`mt-` keeps the
  previous single-colour scatter and says so in the log, instead of a colour bar reading 0% everywhere.
  The gene pattern now lives once in `qcutil.MITO_PREFIXES`/`qcutil.mito_percent()`, shared with
  `filterbycount`, so the number a cell is filtered on and the number colouring it cannot drift apart.
- `.obs` of every matrix from `filterbycount` onwards — including `result/{sample}.h5ad` and
  `result/{sample}_obs.txt.gz` — gains `raw_total_counts`, `raw_n_genes_by_counts` and
  `raw_pct_counts_mt`: the same three QC metrics computed on the uncorrected Cell Ranger counts. The
  unprefixed columns were, and remain, post-ambient-correction; nothing about filtering changes, since no
  threshold is applied to the new columns. `raw_` means *before ambient correction* (source:
  `filtered_feature_bc_matrix.h5`, the same cells), not the all-droplets `raw_feature_bc_matrix.h5`.
  `filterbycount` now takes both matrices as named inputs (`corrected=`, `raw=`) and errors if a corrected
  barcode is missing from the Cell Ranger matrix rather than writing a silently NaN column.

# v0.3.1 - Aug 5, 2026

Packaging release. No behaviour change: the workflow, its outputs and its numbers are identical to v0.3.0.

- The v0.3.0 sdist on PyPI was built one commit before the `v0.3.0` tag, so the published copy of
  `scripts/nuclear_fraction.py` still pointed at `docs/v0.2.0_plan.md 4.2` — a file renamed to
  `docs/design.md` in the tagged tree. PyPI does not allow replacing a released file, so this release
  exists to make the published artifact and the tag agree.
- `envs/cellqc.yaml` floors `r-soupx>=1.6.2`, matching the bioconda recipe, which had carried the pin
  alone. SoupX upstream has been dormant since 1.6.2 (2022-11-01) while conda-forge still carries builds
  back to 1.4.5, so the floor asks for the final release rather than letting a constrained solve pick an
  old one. Verified: the constraint installs 1.6.2 and every SoupX call `ambient.R` makes
  (`autoEstCont(forceAccept=)`, `adjustCounts(roundToInt=)`, `SoupChannel`, `setClusters`) is present.

# v0.3.0 - Aug 5, 2026

## Breaking changes

- Output layout. The final matrix is now `result/{sample}.h5ad` — postproc's output, what a user actually
  takes away — and the pre-integration matrix moved from `result/` to `filterdoublet/{sample}.h5ad`, where
  every other stage's output already lives. The `postproc/` directory is gone. Downstream code that read
  `postproc/{sample}.h5ad` should read `result/{sample}.h5ad`; code that read the old `result/{sample}.h5ad`
  wants `filterdoublet/{sample}.h5ad`, which is now `temp()` — see below. This resolves the open question
  in `docs/design.md` §8.3.
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
  publicly and is purely about testing the package. Gone: `tests/main.sh` (depended on lab-specific
  helper scripts and called `cellqc` without a sample file), `tests/mwe/slurm_cellqc.sh` (a
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

