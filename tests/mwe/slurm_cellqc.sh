#!/usr/bin/env bash
#SBATCH --job-name=cellqc_mwe
#SBATCH --partition=free
#SBATCH --account=ruic20_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.log
#
# End-to-end cellqc v0.2.0 run on the GSE188280 minimal working example.
#
# The `free` partition is preemptible. --requeue plus --open-mode=append lets a
# preempted job restart and append to the same log; Snakemake picks up where it
# left off because completed outputs are already on disk, so a requeue costs
# only the interrupted rule.

set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate cellqc_v0.2.0

INDIR=/dfs3b/ruic20_lab/jinl14/mrrdir/.local/github/cellqc/tests/CellQC_mwe/cellqc_input/GSE188280_GSM5676874_0715_Macula_Retina
WORKDIR="${WORKDIR:-$PWD}"
OUTDIR="${OUTDIR:-$WORKDIR/out}"

cd "$WORKDIR"

printf 'sample\tcellranger\tnreaction\n' > samples.txt
printf 'GSE188280_Macula\t%s\t1\n' "$INDIR" >> samples.txt

cat > config.yaml <<'YAML'
seed: 42
ambient:
  method: soupx
  compare: [decontx]
filterbycount:
  mincount: 500
  minfeature: 300
  mito: 10
doublet:
  run: [doubletfinder, scdblfinder]
  decider: doubletfinder
  findpK: false
  pK: 0.01
  numthreads: 5
  rate: 0.1
  capacity: 13000
nuclear_fraction:
  numthreads: 12
YAML

echo "=== cellqc $(cellqc --version) starting $(date) on $(hostname) (attempt ${SLURM_RESTART_COUNT:-0}) ==="
cellqc -d "$OUTDIR" -t 16 -c config.yaml -- samples.txt
echo "=== finished $(date) ==="
