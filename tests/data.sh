#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug
outdir=$(mrrdir.sh)/test_data
mkdir -p "$outdir" && cd "$outdir" 
wget -O scPred_reference.rds https://bcm.box.com/shared/static/jtmx0lgvsxhcht0ycsccnumx3v2y8ett.rds
wget -O cellqc_test_data.zip https://bcm.box.com/shared/static/19d41zw8tdcf0p82rs90607paklf54ym.zip
unzip cellqc_test_data.zip

cat > samples.txt <<EOF
sample	cellranger
AMD1	"$outdir/cellqc_test_data/AMD1"
AMD2	"$outdir/cellqc_test_data/AMD2"
EOF
AMD1	/tmp/cellqc/cellqc_test_data/AMD1

cat > config.yaml <<EOF
# samples with Cell Ranger output directories
samples: samples.txt

## configuration for dropkick
dropkick:
  method: multiotsu
  numthreads: 1

## configuration for DoubletFinder
doubletfinder:
  findpK: false
  numthreads: 5
  pK: 0.005

## configuration for scPred
scpred:
  reference: "$outdir"/scPred_reference.rds
  threshold: 0.9
EOF
