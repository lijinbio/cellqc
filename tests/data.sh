#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

set -xe
mkdir -p /tmp/cellqc && cd /tmp/cellqc
wget -O scPred_reference.rds https://bcm.box.com/shared/static/jtmx0lgvsxhcht0ycsccnumx3v2y8ett.rds
wget -O cellqc_test_data.zip https://bcm.box.com/shared/static/19d41zw8tdcf0p82rs90607paklf54ym.zip
unzip cellqc_test_data.zip
