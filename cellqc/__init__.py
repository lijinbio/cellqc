"""cellqc -- QC pipeline for single-cell RNA-Seq data.

The version is declared once, in `pyproject.toml`, and read back here from the
installed distribution metadata so `from cellqc import __version__` keeps
working for the CLI's `--version` and for the two reports.

Consequence worth knowing: `importlib.metadata` reads what `pip install` wrote,
not what the working tree says. In an editable install the number is frozen at
install time, so **bumping the version in pyproject.toml requires a
`pip install -e .` before the reports stamp the new one**.
"""

from importlib.metadata import PackageNotFoundError, version as _version

try:
	__version__ = _version('cellqc')
except PackageNotFoundError:  # running from a source tree that was never installed
	__version__ = '0.0.0+unknown'
