# cellqc — design notes

The rationale behind the pipeline: why each step exists, what it deliberately does not do, and what was
measured rather than assumed. Written for v0.2.0 and kept current where later releases superseded it —
supersessions are marked inline rather than edited away, because the reasoning is the point.

The task board and the installation section that used to live here are gone: the first was
project management for a release that shipped, the second duplicates README.md, which is the
source of truth for installation.

`tests/validate_nuclear_fraction.py` and `cellqc/scripts/nuclear_fraction.py` cite section 2 of this file
for the nuclear-fraction definition and its acceptance gate.

---

## 0. Scope and rationale

v0.2.0 is a deliberate simplification and modernization of the pipeline:

1. Drop dropkick and scPred entirely.
2. Make the install a single conda environment (drop SeuratDisk; prefer conda over GitHub).
3. Add a presentation-ready PDF slide report (beamer/tectonic) alongside the HTML report.
4. Add a nuclear-fraction vs. log10(UMI) scatter plot.
5. Move analysis steps to Python/scanpy where equivalent; keep `.h5ad` as the single interchange format.

### Breaking changes (must be in CHANGELOG)

- `dropkick` and `scpred` config sections are removed. Configs containing them will warn and ignore, not fail.
- `.h5seurat` intermediates and `result/*.h5seurat` are gone. `.h5ad` is the only matrix format.
- Seurat is unpinned from v4 (SeuratDisk removal lifts the constraint); v0.2.0 targets Seurat 5.
- Numerical results for the DoubletFinder step will not reproduce v0.1.0 bit-for-bit, because v0.1.0 set no random seed (see §5.3). This is a fix, not a regression, but old and new runs are not directly comparable.

### Environment facts established before planning (verified, not assumed)

| Component | Finding |
|---|---|
| Single conda env (py3.12 + R 4.5 + all deps) | Solves cleanly, 486 packages — verified by `mamba env create --dry-run` |
| `r-soupx`, `bioconductor-dropletutils`, `bioconductor-zellkonverter`, `r-seurat 5.5.1`, `pysam`, `tectonic` | All on conda-forge/bioconda |
| **DoubletFinder** | **Not on conda.** GitHub-only — the one remaining `remotes::install_github` |
| `lijinbio/DoubletFinder` fork | Obsolete. Upstream v2.0.6 already dispatches on Seurat version via `LayerData(seu, layer="counts")`; the fork is *behind* upstream. Switch to upstream. |
| `DropletQC` | GitHub-only. Reimplementing its one used function in pysam removes this dependency (§4). |
| Test sample | Cell Ranger 10.0.0, 3' v3 (polyA), GRCh38, 38,606 features × 2,139,892 raw barcodes, 13,559 called cells |
| BAM tags | `CB` and `RE` (values E/N/I) both present in the CR 10.0.0 BAM — confirmed by `samtools view` on the first 200k records |

### Doublet calling: two callers, DoubletFinder decides

Both DoubletFinder and scDblFinder run on the same post-filter matrix. **DoubletFinder alone determines
which cells are removed**; scDblFinder's score and class are carried in `.obs` and cross-tabulated in the
report as a concordance check. Rationale:

- scDblFinder is on bioconda and actively maintained (Bioconductor 3.22, v1.24.0), so it costs no new
  install and little runtime (it is substantially faster than DoubletFinder on ~13k cells).
- The two methods are the top performers in the Xi & Li (2021, *Cell Systems*) benchmark, but they are not
  redundant: they use different artificial-doublet strategies (DoubletFinder: kNN on a PCA of real +
  simulated cells; scDblFinder: gradient-boosted classifier on cluster-aware simulated doublets). Agreement
  between them is genuine evidence; disagreement is a flag for manual review.
- Keeping the *decision* with one caller avoids an undeclared ensemble rule. A union would raise the removal
  rate above the assumed multiplet rate; an intersection would depress it. Either is defensible but must be
  chosen deliberately, not by accident — so the default stays single-caller and the second is diagnostic.

The concordance table (2×2, DoubletFinder × scDblFinder) and Cohen's κ go in both reports. If κ is low
(< ~0.5) on a given sample that is a signal the doublet call is unstable for that sample, and the report
will say so rather than presenting the DoubletFinder call as settled.

**Built to be swapped.** Since scDblFinder may later replace DoubletFinder outright, the step is
parameterized rather than hard-wired:

```yaml
doublet:
  run: [doubletfinder, scdblfinder]   # which callers to execute
  decider: doubletfinder              # which one's call removes cells
```

(`skip: false` was also in this section as shipped in v0.2.0. It was removed afterwards: leaving a caller
out of `run` already expresses "do not run it", and there is no useful meaning for "detect doublets and
keep them".)

`filterdoublet.py` reads `decider` and applies that caller's `*_metadata.txt.gz`; every caller in `run`
contributes its score/class to `.obs` regardless. Switching to scDblFinder-only later is then
`run: [scdblfinder], decider: scdblfinder` — a config edit, no code change and no rule rewiring. The
`.obs` column names are namespaced (`doubletfinder_pANN`, `doubletfinder_class`, `scdblfinder_score`,
`scdblfinder_class`) so downstream code never has to guess which caller produced a column.

**On deciding which is "superior" (§8.1).** This needs stating up front: **on real data with no ground
truth, neither caller can be shown superior.** Concordance, removal rate, and how clean the UMAP looks are
all consistency measures, not accuracy measures — a caller that removes more cells is not thereby more
correct. A defensible comparison needs doublets whose identity is known independently:

- **Genetic multiplexing** — a pooled multi-donor library where demuxlet/souporcell/vireo flags
  cross-donor droplets. This gives real doublet labels, but only heterotypic *cross-donor* ones, so it
  measures a lower bound on sensitivity.
- **Cell-hashing** — HTO-labelled samples, same caveat.
- **Simulation** — computationally combine real cells to make synthetic doublets with known labels. Full
  control, but simulated doublets are easier than real ones and the simulation strategy resembles what both
  methods already do internally, which biases the comparison.

None of these exist for the current test sample. So v0.2.0 will record the per-cell scores from both
callers to make such an evaluation *possible* later, and the report will present concordance as concordance
— explicitly not as evidence that either method is right. Any future switch should rest on one of the three
designs above, and §8.1 tracks that as an open evaluation rather than something this release settles.

### Ambient RNA correction: SoupX, DecontX, CellBender

Short answer: **SoupX + DecontX is cheap and worth doing in v0.2.0. CellBender is not too much, but it does
not belong in the same release** — it is a different kind of tool wearing the same label.

The three are not interchangeable methods for one slot:

| | Input | What it does | Output barcodes | Cost |
|---|---|---|---|---|
| **SoupX** | filtered + raw + clustering | Estimates a global contamination fraction ρ from marker genes, subtracts it | unchanged | seconds–minutes, CPU |
| **DecontX** | filtered (raw optional) | Bayesian mixture, per-cell contamination proportion, returns decontaminated counts | unchanged | minutes, CPU |
| **CellBender** | **raw** | Jointly does **cell calling** *and* ambient removal with a variational autoencoder | **redefined** | ~1–3 h/sample, wants a GPU |

SoupX and DecontX slot into the existing `soupx` position with no structural change: same input, same output
barcode set, corrected counts. Adding DecontX is roughly 60 lines of R plus a config switch, and
`bioconductor-celda` is on bioconda. Low cost, real value — the two use independent signals (SoupX keys on
marker genes assumed off in most cells; DecontX on a cluster-conditional mixture), so agreement between
their contamination estimates is meaningful evidence, and disagreement is a flag.

CellBender is the expensive one, and the cost is structural, not just compute:

1. It **replaces the Cell Ranger cell call**. Its output has a different barcode set, so the nuclear-fraction
   step (keyed on `filtered_feature_bc_matrix/barcodes.tsv.gz`) and the barcode-rank plot both need
   rethinking. It is not a swap at the ambient slot; it is a swap at the ambient *and* cell-calling slots.
   Note this is the same job dropkick was doing, which v0.2.0 is deleting — so adding CellBender partly
   re-opens what goal 1 closes.
2. It pulls **pytorch** into the environment (~2–3 GB), against goal 2's "one simple env". Mitigation:
   isolate it with a Snakemake per-rule `conda:` directive so the base install stays light.
3. It wants a **GPU** and hours per sample. Feasible here (`free-gpu` / `ruic20_lab_gpu`), but it makes the
   pipeline no longer runnable end-to-end on CPU, which is a real change in what cellqc is.
4. Its `--expected-cells` / `--total-droplets-included` need per-sample judgement; wrong values degrade
   results quietly. That conflicts with the "runs on a new sample without hand-tuning" property (goal 3.2).

**Measured cost of each, in this env (not estimated — solved):**

- `bioconductor-celda` (DecontX): **18 packages, 24 MB, solves into the main env.** Essentially free.
- `cellbender`: **does not solve at all.** The bioconda build pins `python=3.7`, which is EOL and
  irreconcilable with this env's Python 3.12:

  ```
  ├─ cellbender =* * is installable and it requires
  │  └─ python =3.7 *
  └─ pin on python =3.12 * is not installable
  ```

  So CellBender cannot share the environment even in principle — it needs its own env regardless of the
  pytorch question. (PyPI `cellbender` 0.3.2 supports newer Python; the bioconda recipe is stale. Either way
  it is a separate env.)

**Recommendation:** SoupX + DecontX in v0.2.0 behind `ambient.method` — the marginal cost is 24 MB;
CellBender deferred to v0.3.0 with an isolated env and a GPU Slurm profile, where the cell-calling
interaction can be designed properly instead of bolted on.

```yaml
ambient:
  method: soupx          # soupx | decontx | none   (cellbender in v0.3.0)
  compare: [decontx]     # also run these, report their estimates, do NOT apply them
```

**Important guard.** `compare` reports contamination estimates side by side; it does **not** let anyone pick
whichever correction makes the downstream result look best. Only `method` touches the counts, it is recorded
in the report and in `.uns`, and choosing it post hoc by inspecting downstream results would be exactly the
tune-until-it-looks-good failure this pipeline exists to prevent. The comparison is a diagnostic — if SoupX
and DecontX disagree sharply on a sample, that is information about the sample, not a menu.

Task additions: **C6** `scripts/decontx.R` + `ambient.method` dispatch — **DONE**;
**C7** ambient-estimate comparison panel in both reports — **DONE**;
**G1** CellBender, deferred to v0.3.0 — **DEFERRED to v0.3.0**.

### Explicitly rejected alternatives

- **Ensemble/consensus doublet removal (union or intersection of the two callers)** — rejected as the
  default for the reason above. The data to build one is now recorded per cell, so it can be revisited with
  evidence.
- **`sceasy`/reticulate for h5ad→R** — works, but routes R through a Python bridge in the same env. `zellkonverter::readH5AD(reader="R")` is pure rhdf5, no interop, deterministic. Chosen.
- **R writing `.h5ad` back out** — rejected. The R doublet step emits only a small per-barcode metadata TSV; Python applies it to the `.h5ad`. Avoids a full re-write of the matrix through a second serializer and removes a class of layer/dtype-loss bugs.

---

## 1. Target workflow

```
cellranger outs/
  ├─→ ambient (R)            → ambient/{s}.h5 (10x HDF5), contamination estimate + plot
  │     soupx (default) | decontx ; non-applied methods run as `compare` diagnostics
  ├─→ nuclear_fraction (py)  → nuclear_fraction/{s}.txt.gz, NF-vs-log10UMI scatter
  └─→ barcoderank (py)       → barcoderank/{s}.{pdf,png}, knee stats        [new]

ambient/{s}.h5
  └─→ filterbycount (py/scanpy) → filterbycount/{s}.h5ad, violin before/after, filter_ncell.txt

filterbycount/{s}.h5ad
  ├─→ doubletfinder (R)      → doubletfinder/{s}_metadata.txt.gz, doublet_ratio.txt, UMAP/tSNE/pANN plots
  └─→ scdblfinder (R)        → scdblfinder/{s}_metadata.txt.gz, score/class            [new, diagnostic]
      └─→ filterdoublet (py) → filterdoublet/{s}.h5ad  (removal by DoubletFinder; both calls kept in .obs)

filterdoublet/{s}.h5ad + nuclear_fraction/{s}.txt.gz
  └─→ postproc (py)          → result/{s}.h5ad  (the final matrix; see §9b question 3)

all of the above
  ├─→ qcreport (py/Jinja2)   → result/report.html
  └─→ slidereport (py/Jinja2 → LaTeX → tectonic) → result/report_slides.pdf   [new]
```

Rule count goes 11 → 10, but the structural win is elsewhere: the conditional-chaining invariant collapses
from three skip flags to one (`doublet.skip`; removed after v0.2.0 — doublet detection always runs and
the callers are chosen with `doublet.run`, so no skip flag survives). In v0.1.0 the "last enabled stage must emit `result/{s}`"
rule had to be threaded through `filterbycount`, `doubletfinder`, and `scpred` simultaneously; now
`filterdoublet.py` is the single terminus and only one conditional remains. That is the real simplification
from goal 1, and it is why adding DecontX and scDblFinder does not re-complicate the graph — both are
selected *inside* a rule by config, not by adding conditionally-included rule files.

## 2. Nuclear fraction in pysam

### 2.1 Definition (from DropletQC source, `R/nuclear_fraction_tags.R`)

```
nuclear_fraction = n_intronic / (n_intronic + n_exonic)
```

per cell barcode, over reads carrying both `CB` and `RE` tags, restricted to the barcodes in
`filtered_feature_bc_matrix/barcodes.tsv.gz`. Intergenic reads (`RE:A:I`) are excluded from both
numerator and denominator. Exon tag `E`, intron tag `N`.

### 2.2 Implementation and known differences from DropletQC

Single pass over the BAM with `pysam`, parallelized by contig (`multiprocessing`), counting `E` and `N`
per barcode in the cell-barcode set. Deterministic — no seed needed.

Two deliberate differences, both to be quantified, not hidden:

1. **No tile double-counting.** DropletQC splits the genome into 100 tiles and uses `scanBam(which=)`;
   a read spanning a tile boundary is counted in both tiles. A single pass counts each record once.
   Expected effect: negligible (boundary reads are ~1e-5 of the total) but it makes exact equality impossible.
2. **Read filtering.** DropletQC's default `ScanBamParam` applies no flag filter, so secondary and
   duplicate alignments are counted. The reimplementation will match this default so the comparison is
   like-for-like, and additionally record counts under a "primary-only" filter so the sensitivity of the
   statistic to that choice is measurable.

**Acceptance gate (D2).** Compare against the existing DropletQC output at
`tests/CellQC_mwe/cellqc_outdir/.../nuclear_fraction/*.txt.gz` (13,559 barcodes, same sample). Accept only if:

- identical barcode set (13,559);
- Pearson and Spearman r > 0.999;
- median |Δ| < 0.001 and max |Δ| < 0.01.

Report the actual numbers and a Bland–Altman plot of the difference. **If the gate fails, D1 is reverted and
DropletQC stays as a GitHub dependency** — matching the reference implementation matters more than removing
one install. This is recorded as a real possibility, not a formality.

### 2.3 Scatter plot

x = log10(total UMI per called cell, from `filtered_feature_bc_matrix.h5` — i.e. pre-SoupX raw counts, so
the axis means what it says), y = nuclear fraction, one point per called cell. Rasterized scatter at 500 dpi
inside a vector PDF. This is the DropletQC diagnostic view: damaged cells appear as high-NF/low-UMI, empty
droplets as low-NF/low-UMI.

**Interpretation is descriptive only.** v0.2.0 will *not* use NF to filter cells. DropletQC's
`identify_empty_drops`/`identify_damaged_cells` rely on a 2D density model whose thresholds are sample- and
tissue-dependent; applying them automatically across samples without per-sample review would be exactly the
kind of unvalidated auto-filtering that should not be buried inside a pipeline. NF is carried in `.obs` and
plotted so a human can decide.

---

## 3. Analysis-correctness items found while reading v0.1.0

These are pre-existing issues, not new work items invented for the sake of it. Each changes a number.

### 3.1 No random seed anywhere — **fix** — WORSE THAN FIRST ASSESSED

Originally scoped as "the doublet calls are not reproducible". Validation showed the problem reaches the
**count matrix itself**:

- `doubletfinder.R` runs `RunPCA`/`RunUMAP`/`RunTSNE` and `doubletFinder` (which samples artificial
  doublets) with no `set.seed`.
- **`SoupX::adjustCounts(roundToInt=TRUE)` performs *randomised* rounding** —
  `out@x = floor(out@x) + rbinom(length(out@x), 1, out@x - floor(out@x))`, plus `sample()` calls in the
  allocation step. With no seed, **every v0.1.0 run produced a different integer count matrix**, which
  propagates into `pct_counts_mt` and flips cells sitting near the mitochondrial threshold.

This is the confirmed cause of the 1-cell difference between v0.1.0 and v0.2.0 at `filterbycount`
(11,233 vs 11,234) — not a `>=`/`>` boundary bug. 133 cells sit within ±0.1 percentage points of the 10%
mito cutoff, so a ±1 count change is enough to move one across.

Fix: a top-level `seed` key (default 42), set in every stochastic script, and re-seeded **before each
ambient method** so a method's result does not depend on which other methods ran first.

### 3.2 Hard-coded doublet-rate constant — **document, do not change**

`doubletratio = ncol(x) * 0.1 / (nreaction * 13000)` encodes the 10x multiplet-rate rule of thumb
(~0.8% per 1,000 cells recovered). It is reasonable for 3' v3 but is not a law, and it silently assumes a
~13,000-cell target per reaction. v0.2.0 keeps the formula (changing it would break comparability with
existing runs) but exposes the `0.1` and `13000` as named config parameters with the current values as
defaults, and prints the assumed rate in the report so it is visible rather than buried in a script.

### 3.3 Homotypic doublet proportion is never modelled

DoubletFinder's `modelHomotypic()` adjusts `nExp` for doublets formed from same-type cells, which are
undetectable. v0.1.0 uses the unadjusted Poisson `nExp`, which **over-removes** cells. Correcting this
requires cell-type or cluster labels. Since scPred is being removed, the honest options are (a) use Seurat
clusters as the annotation for `modelHomotypic`, or (b) leave it unadjusted and say so. **Recommendation:
(b) for v0.2.0** — leave unadjusted, and state the direction of the bias in the report, because introducing
a clustering step purely to feed the correction adds a tuned parameter (resolution) to a QC pipeline.
Flagged for the user's decision rather than changed unilaterally.

### 3.4 SoupX contamination is estimated but its effect is never reported

`autoEstCont` estimates rho and `adjustCounts` applies it, but nothing quantifies the impact. Add to the
report: rho per sample, and total counts removed as a fraction. Cheap, and it makes an invisible correction
auditable.

---

## 4. PDF slide report

`scripts/slidereport.py` renders `scripts/template/slides.tex.jinja2` and compiles with `tectonic`.

**Template is data-driven, not hard-coded (goal 3.2):** the template loops over the `samples` DataFrame and
over a declared list of figure blocks. Adding a sample requires no template edit; a missing optional figure
(e.g. no BAM → no NF plot) drops its frame rather than failing the build. Jinja2 is configured with LaTeX-safe
delimiters (`((* *))`, `((( )))`) so `{}` in TeX is untouched.

Deck structure, ordered as a narrative rather than a figure dump:

1. Title — run ID, cellqc version, date, seed, sample count
2. Pipeline overview — what ran, what was skipped, and the config actually used
3. **Per sample:**
   1. Cell Ranger metrics table (goal 3.1)
   2. Barcode rank plot with the called-cell count annotated (goal 3.1)
   3. SoupX rho + counts removed
   4. QC violins before/after, with thresholds drawn as lines
   5. Nuclear fraction vs log10 UMI (goal 4)
   6. Doublet calls: pANN violin, UMAP/tSNE
4. Cross-sample summary table — cells at every stage, one row per sample
5. Caveats slide — seed, unmodelled homotypic doublets, NF not used for filtering

Beamer settings per house style: `aspectratio=169`, `10pt`, Madrid theme, `\definecolor{retina}{RGB}{16,110,120}`,
`structure` fg and `frametitle` bg set to it, `\large` frame titles, 12pt text margins, navigation symbols off.

The metrics table has 20 columns and will not fit a 16:9 slide — it will be split into two `tabular`s
(sequencing/mapping | cells/genes) rather than shrunk to unreadable type.

---

## 5. Figure policy (house style)

Every figure: `matplotlib` → PDF with text as font objects (`rcParams['pdf.fonttype']=42`), dense layers
(scatter, hexbin) with `rasterized=True` at `dpi=600`; a PNG twin at 300 dpi for the HTML report only.
(Shipped as 500/200; raised afterwards — 200 dpi was soft under zoom.)
R figures (`ggplot2`) use `ggsave(device=cairo_pdf)` for the same reason. No `useDingbats`.

---

## 6. Validation

Run on `tests/CellQC_mwe/cellqc_input/GSE188280_GSM5676874_0715_Macula_Retina` (CR 10.0.0, 13,559 cells).
v0.1.0 outputs for this sample already exist in `cellqc_outdir/`, so every stage is comparable.

| Check | Criterion |
|---|---|
| Nuclear fraction | §4.2 gate vs DropletQC |
| SoupX | rho within 1e-6 of v0.1.0 (same code path, R only) |
| filterbycount | cell count identical to v0.1.0 (scanpy vs Seurat must agree; any difference is a bug in C1 and must be explained, not accepted) |
| DoubletFinder | cell count within Poisson noise of v0.1.0; exact match not expected (no seed in v0.1.0) |
| scDblFinder | runs, produces a score per cell, removal rate close to the assumed `dbr` |
| Concordance | 2×2 table + κ reported. **No pass/fail threshold** — this is a descriptive statistic, and setting a gate on it would imply an accuracy claim the data cannot support (§ "On deciding which is superior") |
| result `.h5ad` | dims, dtypes, `.obs` columns, no NaN in `pANN`/`nuclear_fraction`, X is integer counts |
| Reports | HTML renders; PDF compiles; text selectable in the PDF |

Slurm: `--partition=free --account=ruic20_lab`, `--requeue` with automatic resubmission on preemption.
The BAM pass is the expensive step (~15 GB, ~233M reads); it gets its own job with the thread count from config.
The submission script is kept in `tests/CellQC_mwe/` for the record.

---

### 6.1 Deferred: DoubletFinder vs scDblFinder evaluation — **DEFERRED**

Blocked on ground-truth data. Requires a genetically multiplexed or cell-hashed library; see the discussion
above. v0.2.0's job is only to record both callers' per-cell scores so this evaluation becomes possible.
Until it is run, "scDblFinder is better" and "DoubletFinder is better" are both unsupported on this data.

---

## 7. Decisions taken (were open questions)

- **§5.3 homotypic doublets — RESOLVED: leave `nExp` unadjusted.** `modelHomotypic()` is not called. The
  consequence must be stated in every report: because same-type doublets are undetectable and no
  annotation-based correction is applied, `nExp` is an *over*-estimate of the detectable doublet count, so
  the step **removes somewhat more cells than the true heterotypic doublet count**. The bias direction is
  known and constant, which is why it is acceptable — it is disclosed, not hidden. Rejected the alternative
  (cluster to feed `modelHomotypic`) because it injects a tuned resolution parameter into QC.
- **§9.2 nuclear fraction — RESOLVED: automatic, no config flag.** Presence of
  `{cellranger}/possorted_genome_bam.bam` (+ `.bai`) is detected per sample at DAG-construction time.
  Present → the step runs; absent → it is skipped for that sample and `postproc` drops the NF input.
  No `skip:` key. Mixed cohorts (some samples with BAM, some without) work, and which samples got NF is
  recorded in both reports so a missing NF column is never silently mistaken for a computed one.

## 8. Remaining open questions

1. **Homotypic doublets.** Leave `nExp` unadjusted (recommended) or add a clustering step to feed
   `modelHomotypic()`? Unadjusted over-estimates the *detectable* doublet count, a bias whose direction is
   known and stated in every report.
2. ~~**Nuclear fraction with no BAM.** v0.1.0 makes it mandatory. v0.2.0 adds `nuclear_fraction.skip` so
   samples without a BAM can run. Confirm that is wanted.~~ **Resolved in v0.2.0:** no skip key at all —
   the step is enabled per sample by the presence of an indexed BAM, detected at DAG-construction time.
3. **`postproc/` vs `result/`** — ~~both are kept as-is. Confirm the split is still useful now that the
   pipeline is shorter.~~ **Resolved after v0.2.0:** the split was not useful. The final matrix is what a
   user wants and it now *is* `result/{s}.h5ad` (postproc's output); the pre-integration matrix moved to
   `filterdoublet/{s}.h5ad`, where every other stage's output already lives. `postproc/` is gone.
