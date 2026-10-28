#' Methods for graphical log-linear fits
#' @param x,object An atom_loglinear_fit.
#' @param type Return compact factors (default) or an explicit probability array.
#' @param force Override full_table_max_cells. R's integer array limit still applies.
#' @param ... Unused.
#' @return fitted returns a factor list or full probability array (summing to
#'   one, not fitted counts). logLik returns the multinomial probability
#'   log-likelihood kernel sum_x n_x log p_hat(x) conditional on the total
#'   sample size, with the model-independent multinomial combinatorial constant
#'   log(n!) - sum_x log(n_x!) omitted and without model df. deviance returns a
#'   numeric saturated-reference deviance, possibly NA. summary returns
#'   diagnostics and engine provenance; print returns x invisibly.
#' @name atom_loglinear_fit
NULL
# Compact per-atom engine provenance table for summary().
.loglinear_engine_provenance <- function(object) {
  fits <- object$local_fits
  field <- function(name, default) vapply(fits, function(z)
    if (is.null(z[[name]])) default else as.character(z[[name]])[[1L]], "")
  reasons <- field("fallback_reason", NA_character_)
  table <- data.frame(atom_id=vapply(fits,`[[`,0L,"atom_id"),
    requested_engine=field("requested_engine", NA_character_),
    attempted_engines=vapply(fits,function(z) paste(z$attempted_engines,collapse=","),""),
    engine=field("engine",""), fallback=vapply(fits,function(z) isTRUE(z$fallback),logical(1)),
    fallback_reason=reasons, stringsAsFactors=FALSE)
  list(requested_engine=unique(field("requested_engine", NA_character_)),
    fallback_count=sum(table$fallback),
    fallback_reasons=unique(reasons[!is.na(reasons)]), table=table)
}
#' @rdname atom_loglinear_fit
#' @method print atom_loglinear_fit
#' @export
print.atom_loglinear_fit <- function(x, ...) {
  cat(sprintf("Graphical log-linear fit: n=%d, variables=%d, atoms=%d\n",
    x$suffstats$n,length(x$levels),length(x$local_fits)))
  cat(sprintf("Cells: full=%g, sum of atoms=%g, largest atom=%g; logLik=%.6f\n",
    x$cells$full,x$cells$total_local,x$cells$max_local,as.numeric(stats::logLik(x))))
  invisible(x)
}
#' @rdname atom_loglinear_fit
#' @method summary atom_loglinear_fit
#' @export
summary.atom_loglinear_fit <- function(object, ...) {
  provenance <- .loglinear_engine_provenance(object)
  c(list(n=object$suffstats$n),decomp_cost(object$decomp),
    list(incomplete_atom_count=sum(!object$decomp$complete),level_counts=lengths(object$levels),
      cells=object$cells,engines=unique(vapply(object$local_fits,`[[`,"","engine")),
      requested_engine=provenance$requested_engine,fallback_count=provenance$fallback_count,
      fallback_reasons=provenance$fallback_reasons,engine_provenance=provenance$table,
      convergence=object$convergence,logLik=stats::logLik(object),
      deviance=stats::deviance(object),timing=object$timing))
}
#' @rdname atom_loglinear_fit
#' @export
logLik.atom_loglinear_fit <- function(object, ...) object$likelihood$logLik
#' @rdname atom_loglinear_fit
#' @export
deviance.atom_loglinear_fit <- function(object, ...) object$likelihood$deviance
#' @rdname atom_loglinear_fit
#' @export
fitted.atom_loglinear_fit <- function(object, type=c("factorized","full"), force=FALSE, ...) {
  type <- match.arg(type)
  if (type=="factorized") object$factorized else .loglinear_full(object,force)
}
