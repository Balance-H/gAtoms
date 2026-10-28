#' Fit a graphical log-linear model from local categorical margins
#' @param data A complete categorical data.frame, or a loglinear_suffstats
#'   created by margin_counts with every ordered atom and separator margin.
#' @param graph A supported simple undirected graph with variable names.
#' @param decomp Optional compatible atom_decomp. It is validated, not recomputed.
#' @param engine Incomplete-atom engine: \code{"gRim"} (default, uses
#'   \code{gRim::effloglin}) or \code{"loglin"} (uses \code{stats::loglin}).
#'   The selected engine is used as given and a failure is reported as an error;
#'   there is no automatic fallback between engines. Complete atoms always use
#'   the closed form. Each local fit records requested_engine, attempted_engines,
#'   the successful engine, fallback, fallback_reason and fallback_message; with
#'   an explicit engine the attempted engine is the requested one.
#' @param control Named list: tol=1e-8 (engine convergence threshold),
#'   contract_tol=1e-6 (margin error divided by n), maxit=1000,
#'   local_table_max_cells=1e6 and full_table_max_cells=1e6.
#' @return An atom_loglinear_fit with graph, decomp, levels, suffstats,
#'   local_fits, separator_fits, factorized probabilities, likelihood,
#'   convergence, diagnostics, cell counts, timing, engine, call and version.
#' @details Variables must exactly match graph names; raw columns are reordered.
#'   Factors retain declared unused levels. Character and integer/logical labels
#'   have deterministic frozen levels. Missing values and continuous numeric
#'   columns are rejected. Sampling zeros are supported; structural-zero
#'   constraints are not implemented.
#'
#'   Only atom and separator count arrays are allocated. The fitted probability
#'   is the product of atom factors divided by separator factors, componentwise.
#'   A zero separator marginal implies zero support; its 0/0 configurations are
#'   assigned probability zero. All fitted generator margins are checked.
#'   Engine iteration counts are NA because these APIs do not return them;
#'   convergence is established from the sufficient-margin equations.
#'
#'   logLik() reports the multinomial probability log-likelihood kernel
#'   sum_x n_x log p_hat(x), conditional on the total sample size. The
#'   multinomial combinatorial constant log(n!) - sum_x log(n_x!) is omitted;
#'   it is independent of the fitted model and therefore cancels in
#'   likelihood-ratio deviances. Saturated likelihood uses the same kernel.
#'   Zero-count cells contribute zero and are absent from the sum. With
#'   sufficient margins lacking observed_counts, deviance is NA with a
#'   diagnostic. Model parameter degrees of freedom and AIC/BIC are not
#'   provided. Local failures raise atom_loglinear_error with atom_id,
#'   vertices, engine and status; the error also carries attempted_engines and
#'   engine_errors for the requested engine. No partial successful fit is
#'   returned.
#' @examples
#' G <- igraph::set_vertex_attr(igraph::make_ring(3, circular=FALSE),
#'                             "name", value=c("A","B","C"))
#' dat <- expand.grid(A=c("a","b"), B=c("x","y"), C=c("u","v"))
#' fit <- fit_loglinear_graph(dat, G)
#' fit
#' stats::logLik(fit)
#' stats::deviance(fit)
#' stats::fitted(fit)
#' @export
fit_loglinear_graph <- function(data, graph, decomp=NULL,
    engine=c("gRim","loglin"), control=list()) {
  start <- proc.time()[[3L]]; call <- match.call()
  engine <- match.arg(engine); control <- .loglinear_control(control)
  norm <- normalize_graph(graph); graph <- norm$graph
  if (!norm$named || !igraph::vcount(graph)) stop("Log-linear fitting requires a nonempty graph with variable names.",call.=FALSE)
  variables <- norm$vertex_map$label
  prepared <- if (inherits(data,"loglinear_suffstats"))
    .validate_loglinear_suffstats(data,variables,control$local_table_max_cells) else
    .validate_loglinear_data(data,variables)
  if (is.null(decomp)) {
    .t0 <- proc.time()[[3L]]
    decomp <- atom_decompose(graph)
    tdec <- unname(proc.time()[[3L]] - .t0)
  } else {
    if (!inherits(decomp,"atom_decomp")) stop("decomp must be an atom_decomp.",call.=FALSE)
    valid <- validate(decomp,expensive=TRUE)
    if (!isTRUE(valid)) stop("Invalid decomposition: ",paste(valid$errors,collapse="; "),call.=FALSE)
    if (!identical(norm$vertex_map,decomp$vertex_map) ||
        !isTRUE(all.equal(as_dgCMatrix(graph),as_dgCMatrix(decomp$graph))))
      stop("decomp must match graph edges, variable names and order.",call.=FALSE)
    tdec <- 0
  }
  t <- proc.time()[[3L]]
  sets <- c(decomp$atoms_index,decomp$separators_index)
  ss <- if (inherits(prepared,"loglinear_suffstats")) prepared else
    .loglinear_margins(prepared,sets,control$local_table_max_cells)
  keys <- vapply(.loglinear_sets(ss$sets,variables),.set_key,"")
  requested <- match(vapply(sets,.set_key,""),keys)
  if (anyNA(requested)) stop("Sufficient statistics lack a required atom or separator margin.",call.=FALSE)
  tmargins <- proc.time()[[3L]]-t; t <- proc.time()[[3L]]
  k <- length(decomp$atoms_index)
  local_fits <- lapply(seq_len(k),function(i) {
    v <- sort(decomp$atoms_index[[i]])
    .fit_loglinear_atom(ss$margins[[requested[i]]],igraph::induced_subgraph(graph,v),
      decomp$complete[i],engine,control,i,variables[v])
  })
  separator_ids <- requested[k+seq_along(decomp$separators_index)]
  cached <- lapply(unique(separator_ids),function(j)
    .fit_loglinear_separator(ss$margins[[j]],ss$sets[[j]],j))
  separator_fits <- cached[match(separator_ids,unique(separator_ids))]
  tlocal <- proc.time()[[3L]]-t; t <- proc.time()[[3L]]
  factors <- .loglinear_factorized(local_fits,separator_fits,decomp,ss$levels,ss$n,control$contract_tol)
  term <- function(z) .loglinear_count_loglik(z$observed,z$fitted_prob)
  ll <- sum(vapply(local_fits,term,0.0))-sum(vapply(separator_fits,term,0.0))
  sat <- dev <- NA_real_; diagnostic <- "Observed configuration counts unavailable; saturated likelihood and deviance unavailable."
  if (!is.null(ss$observed_counts)) {
    sat <- .loglinear_count_loglik(ss$observed_counts,ss$observed_counts/ss$n)
    dev <- 2*(sat-ll); diagnostic <- NULL
    if (abs(dev)<control$contract_tol*max(1,ss$n)) dev <- 0
  }
  cells <- list(full=.loglinear_cells(ss$level_counts),log_full=sum(log(ss$level_counts)),
    total_local=sum(vapply(local_fits,function(z) length(z$observed),0)),
    max_local=max(vapply(local_fits,function(z) length(z$observed),0)))
  structure(list(call=call,graph=graph,decomp=decomp,levels=ss$levels,suffstats=ss,
    local_fits=local_fits,separator_fits=separator_fits,factorized=factors,
    likelihood=list(logLik=structure(ll,class="logLik",nobs=ss$n),saturated=sat,deviance=dev),
    convergence=list(converged=TRUE,status=vapply(local_fits,`[[`,"","status")),
    diagnostics=list(deviance=diagnostic,control=control,
      cell_count=if (is.infinite(cells$full)) "Full cell count overflow; use log_full." else NULL),
    cells=cells,timing=list(decomposition=unname(tdec),margins=unname(tmargins),
      local_fitting=unname(tlocal),composition=unname(proc.time()[[3L]]-t),total=unname(proc.time()[[3L]]-start)),
    engine=engine,version="0.1.0"),class="atom_loglinear_fit")
}
