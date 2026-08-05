import os

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


# Callers that run, and the metadata each one emits. The deciding caller is a
# config value, so switching which method removes cells is a config edit rather
# than a rewiring of the DAG. The schema requires at least one caller: doublet
# detection has no skip flag.
doublet_callers = list(config['doublet']['run'])
ambient_methods = [config['ambient']['method']] + list(config['ambient']['compare'])
ambient_methods = [m for m in ambient_methods if m != 'none']


def doublet_metadata(wildcards):
  return [f"{caller}/{wildcards.sample}_metadata.txt.gz" for caller in doublet_callers]


def postproc_input(wildcards):
  """QC'd matrix plus the nuclear fraction, when this sample has a BAM."""
  ins = {'h5ad': f"filterdoublet/{wildcards.sample}.h5ad"}
  if samples.loc[wildcards.sample, 'has_bam']:
    ins['nf'] = f"nuclear_fraction/{wildcards.sample}.txt.gz"
  return ins


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
    'doublet_summary': expand("result/{sample}_doublet_summary.txt", sample=ids),
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
  """Everything `rule all` should ask for, given which steps apply."""
  ids = samples['sample'].tolist()
  # result/{s}.h5ad is postproc's output, the final matrix; filterdoublet/{s}.h5ad
  # is the same cells before the integration prep. Both carry an _obs/_var dump.
  targets = expand(
    [
      "result/{sample}.h5ad", "result/{sample}_obs.txt.gz", "result/{sample}_var.txt.gz",
      "filterdoublet/{sample}.h5ad",
      "filterdoublet/{sample}_obs.txt.gz", "filterdoublet/{sample}_var.txt.gz",
      ],
    sample=ids)
  targets += expand(["barcoderank/{sample}_barcoderank.pdf"], sample=ids)
  targets += expand(["nuclear_fraction/{sample}.txt.gz"], sample=nf_samples)
  targets += ["result/report.html", "result/report_slides.pdf", "result/metrics.csv"]
  return targets
