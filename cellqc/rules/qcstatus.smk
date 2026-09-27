# Where per-sample failures are resolved into a cohort result.
#
# Every per-sample step writes <rule>/{sample}_status.tsv and never fails the
# run on a sample's account (cellqc/qcstatus.py, scripts/guard.R). This
# checkpoint reads them all and writes
#
#   result/qc_status.csv  every sample x step: ok / fallback / failed / skipped
#   result/manifest.tsv   one row per sample: included in result/ or not, and why
#
# It is a checkpoint because which samples reach result/ is only known once
# their steps have run: final_targets() asks for `publish` outputs only for the
# samples the manifest includes.
checkpoint qcstatus:
  input:
    status=all_status_files(),
    summary=expand("filterdoublet/{sample}_doublet_summary.txt", sample=samples['sample'].tolist()),
  output:
    csv="result/qc_status.csv",
    manifest="result/manifest.tsv",
  params:
    samples=samples['sample'].tolist(),
    steps={s: sample_steps(s) for s in samples['sample'].tolist()},
    has_bam=samples['has_bam'].to_dict(),
    ambient_compare=list(config['ambient']['compare']),
  script:
    "../scripts/qcstatus.py"


# Hard links rather than copies: the final matrix is written once. The two
# doublet statistics files stay in filterdoublet/ too (same inode, no extra
# space); postproc's temp() names are removed by Snakemake after this runs.
rule publish:
  input:
    h5ad="postproc/{sample}.h5ad",
    obs="postproc/{sample}_obs.txt.gz",
    var="postproc/{sample}_var.txt.gz",
    summary="filterdoublet/{sample}_doublet_summary.txt",
    concordance="filterdoublet/{sample}_doublet_concordance.txt",
  output:
    h5ad="result/{sample}.h5ad",
    obs="result/{sample}_obs.txt.gz",
    var="result/{sample}_var.txt.gz",
    summary="result/{sample}_doublet_summary.txt",
    concordance="result/{sample}_doublet_concordance.txt",
  run:
    for key in ("h5ad", "obs", "var", "summary", "concordance"):
      publish_file(input[key], output[key])
