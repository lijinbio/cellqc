# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:
"""Per-sample, per-step status: how one bad sample stops being fatal to a cohort.

Before v0.3.4 a sample that failed any step -- SoupX finding no marker genes in
a low-complexity channel, DoubletFinder handed a handful of cells -- killed the Snakemake
run for every sample. Now every per-sample step runs inside `run()` (and its R
twin, `scripts/guard.R`), which always leaves the step's declared outputs behind
plus a status file `<rule>/<sample>_status.tsv`:

	sample<TAB>step<TAB>status<TAB>message

`status` is one of

	ok        the step ran and its outputs are real
	fallback  the step ran, but on a substitute (e.g. uncorrected counts after
	          SoupX failed); the outputs are real and the message says what
	          was substituted
	failed    the step raised; the message is the error
	skipped   a step it needs was not usable; the message is that step's reason,
	          so the root cause travels down the chain instead of being lost

The first row is the step itself. A step may append rows for parts that can
fail on their own without failing it -- `ambient.decontx` for a comparison-only
ambient estimate.

A failed or skipped step still has to produce every declared output, because
Snakemake cannot express an optional one and a missing output would fail the
run. Those are 0-byte placeholders, never mistaken for data: every consumer
reads the status first, and the reports treat an empty file as not computed. A
matrix a failed step may have half-written is truncated to a placeholder too.
Nothing that stood in for a matrix ever reaches `result/`: the `qcstatus`
checkpoint decides which samples are published there.

Only the step body is guarded. Imports and library loads stay outside, so a
broken environment -- which would fail every sample the same way -- still stops
the run instead of being recorded once per sample as a sample problem.
"""

import os
import sys
import traceback

OK = 'ok'
FALLBACK = 'fallback'
FAILED = 'failed'
SKIPPED = 'skipped'
STATUSES = (OK, FALLBACK, FAILED, SKIPPED)
# A step whose outputs downstream steps may consume.
USABLE = (OK, FALLBACK)

COLUMNS = ('sample', 'step', 'status', 'message')
# Long enough for an R error with its call, short enough for a table cell. The
# full error and traceback are in the Snakemake log.
MAX_MESSAGE = 500

# Outputs that must not survive a failure half-written: a truncated count
# matrix that happens to open is worse than none.
MATRIX_SUFFIXES = ('.h5ad', '.h5')


def clean_message(message):
	"""One line, no tabs, bounded: the status file is a TSV read without quoting."""
	text = ' '.join(str(message).split())
	if len(text) > MAX_MESSAGE:
		text = text[:MAX_MESSAGE - 3] + '...'
	return text


def write_status(path, rows):
	"""Write `rows` (dicts with COLUMNS) as the status TSV."""
	with open(path, 'w', encoding='utf-8') as fh:
		fh.write('\t'.join(COLUMNS) + '\n')
		for r in rows:
			fh.write('\t'.join(clean_message(r.get(c, '')) for c in COLUMNS) + '\n')


def read_status(path):
	"""Rows of a status file as dicts, or None when it is absent or empty."""
	if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
		return None
	with open(path, encoding='utf-8') as fh:
		lines = [ln.rstrip('\n') for ln in fh if ln.strip()]
	header = lines[0].split('\t')
	rows = []
	for ln in lines[1:]:
		vals = ln.split('\t')
		vals += [''] * (len(header) - len(vals))
		rows.append(dict(zip(header, vals)))
	return rows or None


def reason(row):
	"""Why a step is unusable, phrased so it can be passed further downstream.

	A skipped step's message already names the root cause, so it is passed on
	unchanged rather than wrapped again: the fourth step in a chain says
	`filterbycount failed: ...`, not `skipped because skipped because ...`.
	"""
	if row['status'] == SKIPPED:
		return row['message'] or f"{row['step']} skipped"
	msg = f": {row['message']}" if row['message'] else ''
	return f"{row['step']} {row['status']}{msg}"


def blockers(paths):
	"""Reasons the required upstream steps in `paths` are unusable (empty if none)."""
	out = []
	for p in paths:
		rows = read_status(p)
		if rows is None:
			out.append(f'no status recorded in {p}')
		elif rows[0]['status'] not in USABLE:
			out.append(reason(rows[0]))
	return out


def usable(path):
	rows = read_status(path)
	return rows is not None and rows[0]['status'] in USABLE


def as_list(files):
	"""A Snakemake named input as a list, whether it holds one path or several."""
	if files is None:
		return []
	if isinstance(files, str):
		return [files]
	return list(files)


class Step:
	"""What a step body reports about itself, beyond raising or returning."""

	def __init__(self, sample, step):
		self.sample = sample
		self.step = step
		self.status = OK
		self.message = ''
		self.extra = []

	def fallback(self, message):
		"""The step completed on a substitute. Downstream steps still run."""
		self.status = FALLBACK
		self.message = message
		log(self.sample, self.step, FALLBACK, message)

	def note(self, message):
		"""An `ok` with something a reader should know (e.g. a caller missing)."""
		self.message = message

	def substep(self, step, status, message=''):
		"""A part of this step that can fail without failing it."""
		self.extra.append({'sample': self.sample, 'step': step, 'status': status, 'message': message})
		if status != OK:
			log(self.sample, step, status, message)

	def rows(self):
		head = {'sample': self.sample, 'step': self.step, 'status': self.status, 'message': self.message}
		return [head] + self.extra


def log(sample, step, status, message):
	print(f'[cellqc] {status.upper()} sample={sample} step={step}: {clean_message(message)}',
		file=sys.stderr, flush=True)


def placeholders(outputs, status):
	"""Leave every declared output in place so Snakemake sees the job complete.

	Anything already written by a step that failed later is kept -- filterbycount
	writes its per-criterion counts before discovering that no cell survived, and
	those counts are the explanation -- except a matrix, which is truncated.
	"""
	for f in outputs:
		if status == FAILED and os.path.exists(f) and f.endswith(MATRIX_SUFFIXES):
			open(f, 'w').close()
		elif not os.path.exists(f):
			os.makedirs(os.path.dirname(os.path.abspath(f)), exist_ok=True)
			open(f, 'w').close()


def run(snakemake, step, main, requires=()):
	"""Run `main(step)` for one sample, recording the outcome instead of dying.

	`requires` are the status files of the upstream steps whose outputs `main`
	consumes; if any is unusable, `main` is not called and the step is skipped.
	`main` receives a `Step` for reporting a fallback, a note or a substep.

	The status file is the output named `status`; every other output gets a
	placeholder if the step did not write it. Always returns normally, so the
	Snakemake job succeeds and the other samples carry on.
	"""
	sample = snakemake.params['sampleid']
	status_file = snakemake.output['status']
	outputs = [f for f in snakemake.output if f != status_file]
	st = Step(sample, step)

	blocked = blockers(as_list(requires))
	if blocked:
		st.status, st.message = SKIPPED, '; '.join(blocked)
		log(sample, step, SKIPPED, st.message)
	else:
		try:
			main(st)
		except Exception as e:
			# The full traceback goes to the Snakemake log; the status table gets
			# the one line a reader needs to decide whether to look there.
			traceback.print_exc()
			st.status, st.message = FAILED, f'{type(e).__name__}: {e}'
			log(sample, step, FAILED, st.message)

	placeholders(outputs, st.status)
	write_status(status_file, st.rows())
	return st
