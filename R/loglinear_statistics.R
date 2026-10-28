.loglinear_failure <- function(status, message, atom_id = NA_integer_, vertices = character(), engine = NA_character_) {
  stop(structure(list(message=message, call=NULL, status=status, atom_id=atom_id,
    vertices=vertices, engine=engine), class=c("atom_loglinear_error","error","condition")))
}

.validate_loglinear_data <- function(data, variables = names(data), levels = NULL) {
  if (!is.data.frame(data) || !ncol(data) || !nrow(data))
    stop("data must be a nonempty categorical data.frame.", call.=FALSE)
  nm <- names(data)
  if (is.null(nm) || anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm))
    stop("Data variable names must be nonempty and unique.", call.=FALSE)
  if (length(variables) != length(nm) || !setequal(nm,variables))
    stop("Data variables must match graph vertex names exactly.", call.=FALSE)
  data <- data[,variables,drop=FALSE]
  lev <- vector("list",length(variables)); names(lev) <- variables
  codes <- matrix(0L,nrow(data),length(variables),dimnames=list(NULL,variables))
  for (j in seq_along(variables)) {
    x <- data[[j]]
    if (anyNA(x)) stop("Categorical data must not contain missing values.", call.=FALSE)
    if (is.factor(x)) {
      l <- base::levels(x)
      if (anyNA(l) || anyDuplicated(l)) stop("Factor levels must be unique and nonmissing.",call.=FALSE)
    } else if (is.character(x)) l <- sort(unique(x),method="radix") else if (
      is.logical(x) || ((is.integer(x) || is.double(x)) && !is.object(x) &&
        all(is.finite(x)) && all(x==trunc(x)))) l <- as.character(sort(unique(x))) else
      stop("Use factors, characters, logical values or finite integer-valued categorical labels.",call.=FALSE)
    if (!is.null(levels)) l <- levels[[j]]
    if (!length(l)) stop("Each variable requires at least one level.",call.=FALSE)
    code <- match(as.character(x),l)
    if (anyNA(code)) stop("Unknown category outside frozen levels.",call.=FALSE)
    lev[[j]] <- l; codes[,j] <- code
  }
  list(codes=codes,variables=variables,levels=lev,n=nrow(data),level_counts=lengths(lev))
}

.loglinear_cells <- function(counts) {
  value <- prod(as.double(counts))
  if (!is.finite(value)) Inf else value
}
.loglinear_index <- function(codes, dims) {
  stride <- c(1,cumprod(as.double(dims)))[seq_along(dims)]
  as.vector((codes-1) %*% stride)+1
}
.loglinear_sets <- function(sets, variables) {
  if (!is.list(sets) || !length(sets)) stop("sets must be a nonempty list of variable subsets.",call.=FALSE)
  lapply(sets,function(s) {
    if (is.character(s)) s <- match(s,variables)
    if (!is.numeric(s) || !length(s) || anyNA(s) || any(s!=trunc(s)) ||
        any(s<1 | s>length(variables)) || anyDuplicated(s))
      stop("Each margin set must contain distinct known variables.",call.=FALSE)
    sort(as.integer(s))
  })
}
.loglinear_margins <- function(validated, sets, max_cells) {
  ids <- .loglinear_sets(sets,validated$variables)
  keys <- vapply(ids,.set_key,""); keep <- !duplicated(keys); ids <- ids[keep]
  margins <- lapply(ids,function(s) {
    dims <- validated$level_counts[s]; cells <- .loglinear_cells(dims)
    if (!is.finite(cells) || cells>max_cells || cells>.Machine$integer.max)
      stop("Requested local margin exceeds max_cells; reduce the local state space.",call.=FALSE)
    array(tabulate(.loglinear_index(validated$codes[,s,drop=FALSE],dims),nbins=as.integer(cells)),
      dim=as.integer(dims),dimnames=validated$levels[s])
  })
  names(margins) <- names(ids) <- keys[keep]
  # Delimited integer codes encode only observed rows, without a Cartesian product.
  row_keys <- do.call(paste,c(unname(as.list(as.data.frame(validated$codes))),list(sep=":")))
  observed <- tabulate(match(row_keys,unique(row_keys)))
  structure(list(n=validated$n,variables=validated$variables,levels=validated$levels,
    level_counts=validated$level_counts,sets=lapply(ids,function(s) validated$variables[s]),
    margins=margins,observed_counts=observed,
    cells=list(full=.loglinear_cells(validated$level_counts),
      log_full=sum(log(validated$level_counts)),total_local=sum(lengths(margins)),
      max_local=max(lengths(margins)))),class="loglinear_suffstats")
}

#' Aggregate categorical margins directly from observations
#' @param data A complete categorical data.frame.
#' @param sets List of variable names or one-based index subsets.
#' @param max_cells Maximum cells in any requested margin (default one million).
#' @return A loglinear_suffstats containing n, variables, frozen levels,
#'   level_counts, deduplicated sets and margins, sparse observed configuration
#'   counts and cell diagnostics. Sets are ordered by input variable order.
#' @details No full contingency table is constructed. Factors retain unused
#'   levels; character levels are sorted and integer/logical labels ordered.
#'   Continuous numeric columns and missing values are rejected. Zero counts
#'   are sampling zeros; structural-zero constraints are not supported.
#' @export
margin_counts <- function(data, sets, max_cells = 1e6) {
  if (!is.numeric(max_cells) || length(max_cells)!=1L || !is.finite(max_cells) || max_cells<1)
    stop("max_cells must be positive and finite.",call.=FALSE)
  .loglinear_margins(.validate_loglinear_data(data),sets,max_cells)
}

.loglinear_count_loglik <- function(counts, prob) {
  positive <- counts>0
  if (any(!is.finite(prob[positive]) | prob[positive]<=0))
    .loglinear_failure("solver_failed","A positive observed count has zero or invalid fitted probability.")
  sum(counts[positive]*log(prob[positive]))
}

.validate_loglinear_suffstats <- function(x, variables, max_cells) {
  if (!identical(x$variables,variables) || !is.list(x$levels) ||
      !identical(names(x$levels),variables) || any(lengths(x$levels)<1L) ||
      any(vapply(x$levels,function(l) !is.character(l) || anyNA(l) || anyDuplicated(l)>0L,logical(1))) ||
      !identical(x$level_counts,lengths(x$levels)))
    stop("Invalid sufficient-statistics variables or frozen levels; use graph variable order.",call.=FALSE)
  if (length(x$n)!=1L || !is.numeric(x$n) || !is.finite(x$n) || x$n<1 || x$n!=trunc(x$n))
    stop("Invalid sufficient-statistics sample size.",call.=FALSE)
  ids <- .loglinear_sets(x$sets,variables)
  if (!is.list(x$margins) || length(ids)!=length(x$margins)) stop("Invalid sufficient margins.",call.=FALSE)
  for (i in seq_along(ids)) {
    s <- ids[[i]]; a <- x$margins[[i]]
    if (!is.numeric(a) || !identical(as.integer(dim(a)),as.integer(x$level_counts[s])) ||
        !identical(dimnames(a),x$levels[s]) || length(a)>max_cells ||
        any(!is.finite(a)) || any(a<0 | a!=trunc(a)) || sum(a)!=x$n)
      .loglinear_failure("invalid_counts","Invalid sufficient margin counts or level metadata.")
  }
  if (!is.null(x$observed_counts) && (!is.numeric(x$observed_counts) || anyNA(x$observed_counts) ||
      any(!is.finite(x$observed_counts)) || any(x$observed_counts<=0 | x$observed_counts!=trunc(x$observed_counts)) ||
      sum(x$observed_counts)!=x$n)) stop("Invalid observed configuration counts.",call.=FALSE)
  x
}
