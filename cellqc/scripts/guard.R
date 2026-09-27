# vim: set noexpandtab tabstop=2:
#
# R twin of cellqc/qcstatus.py -- read that docstring for the contract. Same
# status file, same four statuses, same placeholder rule, so the Python and R
# steps are indistinguishable to everything downstream.
#
# Sourced by the R step scripts from the path in params$guard. A script keeps
# its library() calls outside the guard (a missing package is an environment
# problem, and should stop the run) and wraps everything after them:
#
#   cellqc_guard('doubletfinder', requires=snakemake@input[['upstream']], {
#     ...
#   })
#
# The body is evaluated in the caller's environment, so a script written at top
# level keeps working unchanged inside the braces. Inside it, cellqc_fallback(),
# cellqc_note() and cellqc_substep() report what raising cannot.

CELLQC_USABLE=c('ok', 'fallback')
CELLQC_MAX_MESSAGE=500

.cellqc=new.env()

cellqc_clean=function(x) {
	x=trimws(gsub('[[:space:]]+', ' ', paste(as.character(x), collapse=' ')))
	if (nchar(x)>CELLQC_MAX_MESSAGE) x=paste0(substr(x, 1, CELLQC_MAX_MESSAGE-3), '...')
	x
}

cellqc_log=function(step, status, message) {
	cat(sprintf('[cellqc] %s sample=%s step=%s: %s\n', toupper(status), .cellqc$sample, step,
		cellqc_clean(message)), file=stderr())
}

cellqc_read_status=function(path) {
	if (!file.exists(path) || file.size(path)==0) return(NULL)
	utils::read.delim(path, quote='', colClasses='character', na.strings=character(0), comment.char='')
}

cellqc_reason=function(row) {
	if (row$status=='skipped') return(if (nzchar(row$message)) row$message else paste(row$step, 'skipped'))
	paste0(row$step, ' ', row$status, if (nzchar(row$message)) paste0(': ', row$message) else '')
}

cellqc_blockers=function(paths) {
	out=character(0)
	for (p in paths) {
		s=cellqc_read_status(p)
		if (is.null(s) || !nrow(s)) out=c(out, paste('no status recorded in', p))
		else if (!s$status[1] %in% CELLQC_USABLE) out=c(out, cellqc_reason(s[1, ]))
	}
	out
}

cellqc_fallback=function(message) {
	.cellqc$status='fallback'
	.cellqc$message=message
	cellqc_log(.cellqc$step, 'fallback', message)
}

cellqc_note=function(message) {
	.cellqc$message=message
}

cellqc_substep=function(step, status, message='') {
	.cellqc$extra[[length(.cellqc$extra)+1]]=data.frame(
		sample=.cellqc$sample, step=step, status=status, message=cellqc_clean(message),
		stringsAsFactors=FALSE)
	if (status!='ok') cellqc_log(step, status, message)
}

cellqc_write_status=function(path) {
	rows=c(list(data.frame(sample=.cellqc$sample, step=.cellqc$step, status=.cellqc$status,
		message=cellqc_clean(.cellqc$message), stringsAsFactors=FALSE)), .cellqc$extra)
	utils::write.table(do.call(rbind, rows), file=path, sep='\t', quote=FALSE, row.names=FALSE,
		col.names=TRUE)
}

cellqc_placeholders=function(outputs, status) {
	for (f in outputs) {
		if (status=='failed' && file.exists(f) && grepl('\\.(h5ad|h5)$', f)) {
			file.create(f)
		} else if (!file.exists(f)) {
			dir.create(dirname(f), recursive=TRUE, showWarnings=FALSE)
			file.create(f)
		}
	}
}

cellqc_guard=function(step, expr, requires=character(0)) {
	body=substitute(expr)
	env=parent.frame()
	status_file=snakemake@output[['status']]
	outputs=setdiff(unlist(snakemake@output), status_file)

	.cellqc$sample=snakemake@params[['sampleid']]
	.cellqc$step=step
	.cellqc$status='ok'
	.cellqc$message=''
	.cellqc$extra=list()

	blocked=cellqc_blockers(unlist(requires))
	if (length(blocked)) {
		.cellqc$status='skipped'
		.cellqc$message=paste(blocked, collapse='; ')
		cellqc_log(step, 'skipped', .cellqc$message)
	} else {
		tryCatch(
			withCallingHandlers(
				eval(body, env),
				# Print the call stack while it still exists; tryCatch unwinds it.
				error=function(e) {
					calls=vapply(sys.calls(), function(x) paste(deparse(x, nlines=1L), collapse=' '), '')
					cat('Traceback (innermost last):\n', paste0('  ', utils::tail(calls, 15), '\n'),
						file=stderr(), sep='')
				}),
			error=function(e) {
				call=conditionCall(e)
				# The function name only: a deparsed call carries its arguments and is
				# cut mid-expression by the message limit. The traceback has the rest.
				fn=if (is.null(call)) '' else paste(deparse(call[[1]], nlines=1L), collapse='')
				where=if (nzchar(fn)) paste0(' (in ', fn, ')') else ''
				.cellqc$status='failed'
				.cellqc$message=paste0(conditionMessage(e), where)
				cellqc_log(step, 'failed', .cellqc$message)
			})
	}

	cellqc_placeholders(outputs, .cellqc$status)
	cellqc_write_status(status_file)
	invisible(.cellqc$status)
}
