#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:
#
# Regenerate docs/workflow.png from docs/workflow.dot. Run it after changing the
# workflow.
#
#   bash docs/make_figures.sh
#
# Keeping the source in the repo is why the v0.1.0 diagram could go stale for a
# whole release: it existed only as a PNG. Do not hand-edit the PNG.
#
# Needs graphviz (dot).

set -euo pipefail

DPI=300
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

command -v dot >/dev/null || { echo "graphviz 'dot' not found" >&2; exit 1; }

echo "==> docs/workflow.png"
dot -Tpng -Gdpi=$DPI "$here/workflow.dot" -o "$here/workflow.png"
echo "done."
