import os
import shutil

from snakemake.exceptions import WorkflowError

# See the note in config.smk: .smk files cannot import the cellqc package,
# because the Snakefile's own directory shadows it on sys.path.


def get_cellranger(wildcards):
  return samples.loc[wildcards.sample, 'cellrangerdir']


def get_rawh5(wildcards):
  return os.path.join(samples.loc[wildcards.sample, 'cellrangerdir'], 'raw_feature_bc_matrix.h5')


def get_filteredh5(wildcards):
  return os.path.join(samples.loc[wildcards.sample, 'cellrangerdir'], 'filtered_feature_bc_matrix.h5')


def get_nreaction(wildcards):
  return int(samples.loc[wildcards.sample, 'nreaction'])


# The R twin of cellqc/qcstatus.py. Passed to the R steps as a path because a
# `script:` is copied into .snakemake/scripts/ before it runs, so it cannot find
# a sibling file on its own.
GUARD_R = os.path.join(workflow.basedir, 'scripts', 'guard.R')


# Callers that run, and the metadata each one emits. The deciding caller is a
# config value, so switching which method removes cells is a config edit rather
# than a rewiring of the DAG. The schema requires at least one caller: doublet
# detection has no skip flag.
doublet_callers = list(config['doublet']['run'])
ambient_methods = [config['ambient']['method']] + list(config['ambient']['compare'])
ambient_methods = [m for m in ambient_methods if m != 'none']


def doublet_metadata(wildcards):
  return [f"{caller}/{wildcards.sample}_metadata.txt.gz" for caller in doublet_callers]


def doublet_status(wildcards):
  return [f"{caller}/{wildcards.sample}_status.tsv" for caller in doublet_callers]


def postproc_input(wildcards):
  """QC'd matrix plus the nuclear fraction, when this sample has a BAM.

  The nuclear fraction is optional even then: its status comes along so a
  failed nuclear-fraction step means a matrix without it, not a lost sample.
  """
  s = wildcards.sample
  ins = {'h5ad': f"filterdoublet/{s}.h5ad", 'upstream': f"filterdoublet/{s}_status.tsv"}
  if samples.loc[s, 'has_bam']:
    ins['nf'] = f"nuclear_fraction/{s}.txt.gz"
    ins['nf_status'] = f"nuclear_fraction/{s}_status.tsv"
  return ins


def sample_steps(sample):
  """The per-sample rules that run for `sample`, in pipeline order.

  Each writes <rule>/{sample}_status.tsv. nuclear_fraction is absent for a
  sample without a BAM; the qcstatus checkpoint records it as skipped.
  """
  steps = ['barcoderank', 'ambient']
  if samples.loc[sample, 'has_bam']:
    steps.append('nuclear_fraction')
  steps += ['filterbycount'] + doublet_callers + ['filterdoublet', 'postproc']
  return steps


def all_status_files():
  return [f"{step}/{s}_status.tsv" for s in samples['sample'].tolist() for step in sample_steps(s)]


def included_samples():
  """Samples the qcstatus checkpoint put in result/. Only valid after it ran."""
  manifest = pd.read_csv(checkpoints.qcstatus.get().output['manifest'], sep='\t')
  return manifest.loc[manifest['included'].astype(str) == 'True', 'sample'].astype(str).tolist()


def publish_file(src, dst):
  """Hard-link `src` to `dst`, or copy it where a link is impossible.

  copy2 keeps the modification time, as a link does, so neither makes an input
  of the qcstatus checkpoint look newer than its output.
  """
  if os.path.lexists(dst):
    os.remove(dst)
  try:
    os.link(src, dst)
  except OSError:
    shutil.copy2(src, dst)


def report_inputs():
  """Everything both reports read. Shared so the HTML and the slides can never
  drift into describing different runs."""
  ids = samples['sample'].tolist()
  ins = {
    'ambient_contamination': expand("ambient/{sample}_contamination.txt", sample=ids),
    'ambient_plot': expand("ambient/{sample}_ambient.png", sample=ids),
    'barcoderank_plot': expand("barcoderank/{sample}_barcoderank.png", sample=ids),
    'barcoderank_knee': expand("barcoderank/{sample}_knee.txt", sample=ids),
    'filter_ncell': expand("filterbycount/{sample}_filter_ncell.txt", sample=ids),
    'violin_before': expand("filterbycount/{sample}_violin_before.png", sample=ids),
    'violin_after': expand("filterbycount/{sample}_violin_after.png", sample=ids),
    'doublet_summary': expand("filterdoublet/{sample}_doublet_summary.txt", sample=ids),
    'doublet_concordance': expand("filterdoublet/{sample}_doublet_concordance.txt", sample=ids),
    'qc_status': "result/qc_status.csv",
    'manifest': "result/manifest.tsv",
    'nf_table': expand("nuclear_fraction/{sample}.txt.gz", sample=nf_samples),
    'nf_plot': expand("nuclear_fraction/{sample}_nf_umi.png", sample=nf_samples),
  }
  for caller in doublet_callers:
    ins[f'{caller}_ratio'] = expand(f"{caller}/{{sample}}_doublet_ratio.txt", sample=ids)
  if 'doubletfinder' in doublet_callers:
    ins['doubletfinder_pANN'] = expand("doubletfinder/{sample}_pANN.png", sample=ids)
    ins['doubletfinder_umap'] = expand("doubletfinder/{sample}_umap.png", sample=ids)
  if 'scdblfinder' in doublet_callers:
    ins['scdblfinder_score'] = expand("scdblfinder/{sample}_score.png", sample=ids)
  return ins


def final_targets():
  """Everything `rule all` asks for regardless of which samples succeed.

  The reports, metrics and status cover every sample, included or not. The
  final matrices are `published_targets`, which needs the qcstatus checkpoint.
  """
  ids = samples['sample'].tolist()
  targets = expand(["barcoderank/{sample}_barcoderank.pdf"], sample=ids)
  targets += expand(["nuclear_fraction/{sample}.txt.gz"], sample=nf_samples)
  targets += ["result/report.html", "result/report_slides.pdf", "result/metrics.csv",
    "result/qc_status.csv", "result/manifest.tsv"]
  return targets


def published_targets(wildcards):
  """The final matrix of every sample the qcstatus checkpoint included.

  An input function, not a list: which samples reach result/ is known only once
  their steps have run, so Snakemake re-evaluates this after the checkpoint.
  result/{s}.h5ad is postproc's output, published, with its .obs and .var
  alongside. filterdoublet's matrix is temp() and deliberately not asked for --
  requesting it here would keep Snakemake from ever deleting it.
  """
  return expand(
    ["result/{sample}.h5ad", "result/{sample}_obs.txt.gz", "result/{sample}_var.txt.gz",
      "result/{sample}_doublet_summary.txt", "result/{sample}_doublet_concordance.txt"],
    sample=included_samples())


def run_outcome():
  """Log which samples were excluded, and fail the run if none survived.

  Called from `onsuccess`. A sample failing a step is not a run failure -- it is
  recorded and the rest carry on -- but a run that leaves result/ with no matrix
  at all has not done its job, and must not exit 0 as if it had.
  """
  path = "result/manifest.tsv"
  if not os.path.exists(path):
    return
  manifest = pd.read_csv(path, sep='\t', keep_default_na=False)
  inc = manifest[manifest['included'].astype(str) == 'True']
  exc = manifest[manifest['included'].astype(str) != 'True']
  print(f"cellqc: {len(inc)}/{len(manifest)} sample(s) in result/; per-step status in "
    f"result/qc_status.csv", file=sys.stderr)
  for _, r in exc.iterrows():
    print(f"cellqc: EXCLUDED {r['sample']}: {r['reason']}", file=sys.stderr)
  if len(inc) == 0:
    raise WorkflowError(
      f"every sample ({len(manifest)}) failed a required step; no matrix was written to "
      "result/. The reports and result/qc_status.csv say why.")
