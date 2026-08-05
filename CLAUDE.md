# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Version

The tree is **v0.3.0**. `docs/v0.2.0_plan.md` carries the design rationale, the task board, and the
validation results — read it before changing the workflow.

Environment: `envs/cellqc.yaml` plus one GitHub build (DoubletFinder is not on conda).

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate cellqc_v0.2.0
```

### Hard-won gotchas — do not rediscover these

- **Never `import cellqc` inside a `.smk` file.** Snakemake puts the Snakefile's directory
  (`cellqc/`) on `sys.path`. `cellqc/cli.py` used to be `cellqc/cellqc.py` and shadowed the package;
  it was renamed for this reason. `.smk` files keep their own copies of the few helpers they need;
  `script:` steps run in a separate process and can import the package normally.
- **`zellkonverter::readH5AD(reader="R")` names the assay `X`.** scDblFinder's `.checkSCE()` and DecontX
  both hard-require one named `counts`. Rename it in every R script that reads an `.h5ad`.
- **`doubletFinder(reuse.pANN=)` must be `NULL`, not `FALSE`.** Upstream v2.0.6 changed the check to
  `!is.null()`, so `FALSE` takes the reuse branch and dies with "cannot xtfrm data frames".
- **`SoupX::adjustCounts(roundToInt=TRUE)` is stochastic** (`rbinom` on the fractional part). Without a
  seed the corrected count matrix differs every run. Always `set.seed()` before it — this is why v0.1.0
  results were not reproducible at all, not just in the doublet step.
- **Slurm jobs must not use the session scratchpad**: `/tmp` is node-local, so a compute node cannot see
  it. Use a `/dfs3b` path.
- `bioconda`'s `cellbender` pins `python=3.7` and cannot share this environment.

### Validation status (reference sample GSE188280, CR 10.0.0, 13,559 cells)

- Nuclear fraction vs DropletQC: Pearson r = 1.000000, max |Δ| = 0.00057. Gate passed.
- End-to-end: 13,559 → 11,234 (filter) → 10,223 cells. v0.1.0 gave 10,222; the 1-cell difference is the
  unseeded SoupX rounding, confirmed by a seed-repeat experiment.
- DoubletFinder 1,011 (9.00%) vs scDblFinder 1,153 (10.26%), Cohen's κ = 0.759.
- Reproducibility proven: same seed → identical matrix; seed 42 vs 7 → 1,329 counts of 91.3M differ.

## What this is

`cellqc` is a Snakemake pipeline for QC of single-cell RNA-Seq data, distributed as a pip/conda package with a single `cellqc` console entry point. Input is Cell Ranger output directories; output is cleaned count matrices (`.h5ad`) plus an HTML QC report and a PDF slide deck.

## Commands

```bash
pip install -e .                       # dev install (entry point: cellqc=cellqc.cli:main)
cellqc -h
cellqc -d "$outdir" -t 8 -n -- samples.txt          # dry run; also writes outdir/config_<timestamp>.yaml
cellqc -d "$outdir" -t 8 -c config.yaml -- samples.txt
cellqc -d "$outdir" -c config.yaml $(basharr2cmdopts.sh -o -D -- "${define[@]}")  # -D sample:=:X cellranger:=:/path
```

There is no unit-test framework or CI. `tests/` is two scripts, kept minimal and free of lab-specific
paths because it ships publicly (see `tests/README.md`):

```bash
bash tests/dryrun.sh                                     # seconds, no data: DAG builds, outputs as promised
python tests/validate_nuclear_fraction.py <new.txt.gz> <dropletqc_ref.txt.gz> [outdir]
```

`tests/dryrun.sh` builds stub Cell Ranger directories (`--dry-run` only needs the paths to exist). Run it
after touching `config.smk`, the schema, `Snakefile`, or any rule's outputs — and extend it when a rule's
outputs change, since that is what it asserts.

`tests/validate_nuclear_fraction.py` is a hard gate on the pysam nuclear-fraction reimplementation against
DropletQC; if it fails, the correct response is to revert to DropletQC, not to loosen the thresholds.

The reference run (GSE188280 GSM5676874, 13,559 cells) is documented in `tests/README.md`; the lab Slurm
submission that produced it is not in the repo. Its outputs live in
`/dfs3b/ruic20_lab/jinl14/mrrdir/.local/github/cellqc_v020_final/out`, and the v0.1.0 comparison run in
`.../cellqc/tests/CellQC_mwe/cellqc_outdir/`.

`docs/workflow.png` is generated: `bash docs/make_figures.sh` renders it from `docs/workflow.dot`. Do not
hand-edit the PNG — the v0.1.0 diagram went stale for a whole release because it existed only as a PNG.

The analysis stack lives in `envs/cellqc.yaml`, not `pyproject.toml` (pip cannot install R packages).
`pyproject.toml` declares only what the CLI itself imports. The bioconda recipe is maintained separately at
`../bioconda-recipes/recipes/cellqc/meta.yaml`, and it is the only place the R dependencies are declared as
installable — with the exception of DoubletFinder, which is GitHub-only and cannot be a conda dependency at
all (bioconda forbids network access at install time; the fix would be a separate `r-doubletfinder` recipe).

## Architecture

**Two layers.** `cellqc/cli.py` is a thin click CLI that creates the outdir, `chdir`s into it, and shells out to `snakemake --snakefile cellqc/Snakefile --config samplefile=... outdir=... configfile=... nowtimestr=...`. All real logic is in the workflow. Note the CLI passes the config file *both* via `--configfile` and as a `configfile=` config key (the latter so `config.smk` can resolve paths relative to it).

**Config resolution** happens entirely in `cellqc/rules/config.smk`: it merges a hard-coded `default_params` dict under any user-supplied YAML (per-section, per-key), validates against `schemas/config.schema.yaml`, reads the sample TSV into a module-level `samples` DataFrame indexed by `sample`, and dumps the effective config to `config_<timestamp>.yaml` in the outdir. **The defaults in `config.smk` and the documented defaults in README.md are duplicated — change both together.**

**Relative-path convention:** `cellranger` paths in the sample file are resolved relative to the *sample file's* directory (`sampledir`, used by `get_cellranger`/`get_rawh5` in `common.smk`); (v0.1.0 also resolved `scpred.reference` against the *config file's* directory; that section is gone.)

**Stage selection is by config, not by conditional rule inclusion.** v0.1.0 threaded a
"whichever stage is last must emit `result/{sample}.h5seurat`" invariant through three rule files with
inline conditionals. v0.2.0 has a single terminus: `filterdoublet.py` always writes the QC'd matrix and
`postproc.py` always writes the final one. Alternative methods (SoupX vs DecontX, DoubletFinder vs
scDblFinder) are chosen *inside* a rule from config, so adding a method does not change the DAG shape.
The pipeline is:

```
cellranger outs ─┬─→ ambient (R: soupx|decontx) ─→ filterbycount (py) ─┬─→ doubletfinder (R) ─┐
                 ├─→ barcoderank (py)                                  └─→ scdblfinder  (R) ─┤
                 └─→ nuclear_fraction (py, only if BAM present)                               │
                     filterdoublet (py) ←──────────────────────────────────────────────────────┘
                       └─→ filterdoublet/{s}.h5ad ─→ postproc (py) ─→ result/{s}.h5ad
                     qcreport (py) → report.html    slidereport (py→LaTeX→tectonic) → report_slides.pdf
```

The R doublet steps write only a per-barcode metadata TSV; Python applies it. R never rewrites the matrix,
which keeps a second serializer out of the count path.

**`result/` is what a user takes away**, and the final matrix is `postproc`'s output — not an intermediate
in a `postproc/` directory that readers mistook for a leftover. Every stage writes to `<rulename>/`; only
the final matrix, the doublet statistics and the two reports live in `result/`. `result/{s}.h5ad` carries
`_obs.txt.gz`/`_var.txt.gz` dumps written by `qcutil.write_obs_var`, indexed by `barcode`/`gene`, so `.obs`
and `.var` are readable without anndata. Rules use **named** outputs (`h5ad=`, `summary=`, …) rather than
positional indices.

`filterdoublet/{s}.h5ad` is `temp()`. It is the same cells and counts as `result/{s}.h5ad`, differing only
in the barcode prefix, unique var names and the nuclear-fraction columns, so keeping it wrote every count
matrix to disk twice (72 MB per sample on the reference cohort). **Do not add it to `final_targets()`** —
requesting it would stop Snakemake ever deleting it.

**No step has a skip flag.** Doublet detection is selected by `doublet.run` (schema: at least one
caller), so `doublet.skip` is gone — `config.smk` errors on `skip: true` rather than quietly writing a
matrix with doublets removed, and warns on `skip: false`. `Snakefile` includes both caller rule files
unconditionally.

**The nuclear-fraction step has no skip flag either.** `config.smk` detects an indexed
`possorted_genome_bam.bam` per sample at DAG-construction time and populates `nf_samples`; `final_targets()`
only requests it for those, and `postproc_input()` drops the input for the rest.

**Both reports read `cellqc/reportdata.py` and nothing else**, so the HTML and the slides cannot drift into
describing different runs. Figure paths carry an `{ext}` placeholder — HTML takes `png`, slides take `pdf`.

**Rule/script split.** Every `rules/*.smk` is a single rule that declares I/O and `params` and delegates to `scripts/<name>.{R,py}` via Snakemake's `script:` directive. These scripts are not standalone — they read `snakemake@input[[1]]` / `snakemake.params.sampleid` etc. and cannot be run outside the workflow. Outputs are written to `<outdir>/<rulename>/`, and each stage that feeds the report also writes a small stats file (`*_filter_ncell.txt`, `*_doublet_ratio.txt`, `*_contamination.txt`, `*_knee.txt`, `*_doublet_concordance.txt`) that `cellqc/reportdata.py` collects.

**Reports.** `scripts/qcreport.py` renders `template/index.html.jinja2` into a self-contained
`result/report.html`, base64-inlining every PNG and the CSS (the `data_uri*` helpers come from Snakemake's
own report module). `scripts/slidereport.py` renders `template/slides.tex.jinja2` and compiles it with
`tectonic`; that template uses LaTeX-safe Jinja delimiters (`((* *))`, `((( )))`) so TeX braces survive.

**Figure policy** lives in `cellqc/qcutil.py`: every plot is written as a vector PDF with editable text
(`pdf.fonttype=42`) for the slides *and* a 300 dpi PNG for the HTML, with dense scatter layers
`rasterized=True` at `RASTER_DPI` (600). The R scripts hard-code the same PNG dpi in `ggsave`/`png` —
change them together.

## Conventions

- Packaging is PEP 621 `pyproject.toml` (setuptools backend); there is no `setup.py`. New files under
  `cellqc/` must be covered by `[tool.setuptools.package-data]` or they will be missing from the wheel —
  the workflow is data, not modules, so Snakemake reads it off disk at runtime. `MANIFEST.in` only adds the
  extra files an sdist needs (`README.md`, `CHANGELOG.md`, `LICENSE`, `envs/`).
- Version lives only in `cellqc/__init__.py`; `pyproject.toml` reads it via `dynamic = ["version"]`. Bump it
  together with a `CHANGELOG.md` entry.
- Indentation is inconsistent by design of the original author and enforced by vim modelines: `rules/*.smk` use 2 spaces; `cellqc.py`, `config.smk`, and `scripts/*.py` use **tabs** with `tabstop=2`. Match the file you are editing.
