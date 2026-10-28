.fit_loglinear_separator <- function(observed, vertices, separator_id) {
  list(separator_id=separator_id,vertices=vertices,observed=observed,
       fitted_counts=observed,fitted_prob=observed/sum(observed),engine="closed_form",
       method="saturated",iterations=0L,converged=TRUE,status="closed_form")
}

.loglinear_factorized <- function(local_fits, separator_fits, decomp, levels, n, tol) {
  for (j in seq_along(decomp$separator_atoms)) {
    child <- decomp$separator_atoms[j]; sep <- separator_fits[[j]]
    for (i in c(child,decomp$parent[child])) {
      a <- local_fits[[i]]
      s <- match(sep$vertices,a$vertices)
      if (max(abs(base::margin.table(a$fitted_prob,s)-sep$fitted_prob))>tol)
        .loglinear_failure("invalid_generator_margins","Inconsistent canonical separator margin.",i,a$vertices,a$engine)
    }
  }
  factor <- function(z) list(vertices=z$vertices,prob=z$fitted_prob)
  list(variables=names(levels),levels=levels,n=n,
    atoms=lapply(local_fits,factor),separators=lapply(separator_fits,factor),
    component_orderings=decomp$component_orderings,components=decomp$components)
}

.loglinear_logprob <- function(factorized, data=NULL, codes=NULL) {
  if (is.null(codes)) codes <- .validate_loglinear_data(data,factorized$variables,factorized$levels)$codes
  out <- numeric(nrow(codes)); zero <- rep(FALSE,nrow(codes))
  factors <- c(factorized$atoms,factorized$separators)
  signs <- c(rep(1,length(factorized$atoms)),rep(-1,length(factorized$separators)))
  for (i in seq_along(factors)) {
    f <- factors[[i]]; s <- match(f$vertices,factorized$variables)
    prob <- f$prob[.loglinear_index(codes[,s,drop=FALSE],dim(f$prob))]
    zero <- zero | prob==0
    positive <- prob>0
    out[positive] <- out[positive]+signs[i]*log(prob[positive])
  }
  # A zero separator marginal implies zero global support, including 0/0 cells.
  out[zero] <- -Inf
  out
}

.loglinear_full <- function(object, force) {
  cells <- object$cells$full
  if (!is.logical(force) || length(force)!=1L || is.na(force)) stop("force must be TRUE or FALSE.",call.=FALSE)
  if (!is.finite(cells) || cells>.Machine$integer.max)
    stop("Full table exceeds the supported array size, even with force=TRUE.",call.=FALSE)
  if (!force && cells>object$diagnostics$control$full_table_max_cells)
    stop("Full table exceeds full_table_max_cells; use factorized output or explicitly force=TRUE.",call.=FALSE)
  dims <- lengths(object$levels)
  codes <- arrayInd(seq_len(cells),.dim=as.integer(dims))
  array(exp(.loglinear_logprob(object$factorized,codes=codes)),dim=as.integer(dims),dimnames=object$levels)
}
