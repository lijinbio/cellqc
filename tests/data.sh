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
