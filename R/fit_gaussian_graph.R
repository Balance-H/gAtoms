#' Fit a fixed-graph Gaussian concentration model
#' @param data Numeric observations in rows and variables in columns. If both
#'   graph and data are named, columns are matched to graph vertex names;
#'   otherwise columns are positional. Output always follows graph order.
#' @param graph A supported simple undirected graph.
#' @param decomp Optional validated atom_decomp for exactly the same ordered graph.
#' @param engine Incomplete-atom engine: \code{"gRim"} (default) or \code{"ggm"}.
#'   Complete atoms always
#'   use Cholesky closed form, including when the selected package is absent.
#' @param covariance Recursive completion, no full covariance, or a debug inverse.
#' @param mean Estimate the mean or fix it at zero.
#' @param control Named list: tol (1e-9), contract_tol (1e-6), maxit (10000),
#'   method ("ncd", "covips", "conips"), fallback (TRUE). The last three options
#'   apply to gRim only; fallback retries nonconverged NCD with covips. Tolerances
#'   apply after marginal variance scaling. ggm exposes no iteration limit.
#' @return An atom_gaussian_fit with precision (sparse symmetric Matrix),
#'   covariance (matrix or NULL), mean, h, graph, decomp, suffstats, local_fits,
#'   separator_fits, logdet, likelihood, convergence, diagnostics, timing,
#'   engine, call, and version.
#' @details All fits use SSD/n. Engines receive the corresponding correlation
#'   matrix; results are rescaled without a sample-size adjustment. Ordered
#'   atom precision blocks are added and separator blocks subtracted. Recursive
#'   covariance uses canonical separator covariances, never averaged estimates.
#'   Disconnected components have zero cross-covariance and precision.
#'   The Gaussian module currently requires complete finite observations.
#'   gRim supplies both K and Sigma. ggm supplies Shat as Sigma and its adapter
#'   computes K by one Cholesky inversion. Neither covariance is reconstructed
#'   after numerical projection of precision nonedges to zero.
#'
#'   Required local failures raise an atom_gaussian_error condition with status,
#'   atom_id and vertices, before global assembly. Status distinguishes
#'   singular_statistics, nonconverged and solver_failed. Nonconvergence does
#'   not establish mathematical MLE nonexistence. No ridge is added. ggm requires
#'   positive definite local empirical covariance; gRim may handle singular
#'   empirical covariance when the constrained MLE exists.
#'
#'   Parameter counts are p+m plus p when estimating the mean. Deviance uses a
#'   saturated covariance reference with choose(p,2)-m degrees of freedom; it is
#'   NA with a diagnostic if the saturated empirical covariance is singular.
#' @examples
#' set.seed(1)
#' G <- igraph::make_ring(4, circular = FALSE)
#' X <- matrix(rnorm(400), 100, 4)
#' dec <- atom_decompose(G)
#' fit <- fit_gaussian_graph(X, G, decomp = dec)
#' fit
#' summary(fit)
#' stats::logLik(fit)
#' stats::AIC(fit)
#' stats::BIC(fit)
#' precision(fit)
#' stats::vcov(fit)
#' @export
fit_gaussian_graph <- function(data, graph, decomp = NULL,
    engine = c("gRim", "ggm"), covariance = c("recursive", "none", "solve"),
    mean = c("estimate", "zero"), control = list()) {
  start <- proc.time()[[3L]]; call <- match.call()
  engine <- match.arg(engine); covariance <- match.arg(covariance)
  mean <- match.arg(mean)
  control <- .gaussian_control(control)
  norm <- normalize_graph(graph); graph <- norm$graph
  p <- igraph::vcount(graph); m <- igraph::ecount(graph)
  if (!p) stop("Gaussian fitting requires at least one vertex.", call. = FALSE)
  if (is.null(decomp)) {
    .t0 <- proc.time()[[3L]]
    decomp <- atom_decompose(graph)
    tdec <- unname(proc.time()[[3L]] - .t0)
  } else {
    if (!inherits(decomp, "atom_decomp")) stop("decomp must be an atom_decomp.", call. = FALSE)
    valid <- validate(decomp, expensive = TRUE)
    if (!isTRUE(valid)) stop("Invalid decomp: ", paste(valid$errors, collapse = "; "), call. = FALSE)
    other <- normalize_graph(decomp$graph)
    if (!identical(norm$vertex_map, other$vertex_map) || !identical(norm$named, other$named) ||
        !isTRUE(all.equal(as_dgCMatrix(graph), as_dgCMatrix(other$graph))))
      stop("decomp must correspond to the same graph, vertex labels and order.", call. = FALSE)
    tdec <- 0
  }
  t <- proc.time()[[3L]]
  if (!(is.matrix(data) || is.data.frame(data)) || ncol(data) != p)
    stop("data must have one column per graph vertex.", call. = FALSE)
  labels <- norm$vertex_map$label; cn <- colnames(data)
  if (!is.null(cn) && (anyNA(cn) || any(!nzchar(cn)) || anyDuplicated(cn)))
    stop("Data column names must be unique and nonmissing.", call. = FALSE)
  if (norm$named && !is.null(cn)) {
    if (!setequal(cn, labels)) stop("Data column names must match graph vertices.", call. = FALSE)
    data <- data[, match(labels, cn), drop = FALSE]
  }
  colnames(data) <- labels
  ss <- gaussian_suffstats(data, mean = mean)
  tstats <- proc.time()[[3L]] - t
  t <- proc.time()[[3L]]
  locals <- lapply(seq_along(decomp$atoms_index), function(i) {
    v <- sort(decomp$atoms_index[[i]])
    if ((decomp$complete[i] || engine == "ggm") && ss$n - (mean == "estimate") < length(v))
      .gaussian_failure("singular_statistics",
        "Too few observations for the required positive definite local statistics.", i, v)
    .fit_gaussian_atom(ss$S_ml[v, v, drop = FALSE],
      igraph::induced_subgraph(graph, v), ss$n, decomp$complete[i], engine, control, i, v)
  })
  seps <- lapply(seq_along(decomp$separators_index), function(i) {
    v <- decomp$separators_index[[i]]
    .fit_gaussian_separator(ss$S_ml[v, v, drop = FALSE], v, i)
  })
  tlocal <- proc.time()[[3L]] - t
  t <- proc.time()[[3L]]
  .check_gaussian_separators(decomp, locals, seps, control$contract_tol)
  K <- .assemble_gaussian_precision(locals, seps, p, labels)
  tryCatch(Matrix::Cholesky(K, LDL = FALSE), error = function(e)
    .gaussian_failure("solver_failed", "Assembled precision is not positive definite."))
  Sigma <- .assemble_gaussian_covariance(decomp, locals, seps, K, covariance)
  ld <- .gaussian_logdet(locals, seps); trace <- .gaussian_trace(ss, locals, seps)
  ll <- .gaussian_loglik(ss$n, p, ld$precision, trace)
  df <- p + m + if (mean == "estimate") p else 0L
  ll <- structure(ll, class = "logLik", df = df, nobs = ss$n)
  saturated <- if (ss$n - (mean == "estimate") < p) NULL else
    tryCatch(.gaussian_chol(ss$S_ml), atom_gaussian_error = function(e) NULL)
  dev_message <- NULL
  if (is.null(saturated)) {
    dev <- NA_real_; dev_message <- "Saturated empirical covariance is singular; reference deviance unavailable."
  } else {
    sat_ll <- .gaussian_loglik(ss$n, p, -2*sum(log(diag(saturated))), p)
    dev <- 2*(sat_ll - as.numeric(ll))
    if (abs(dev) < control$contract_tol) dev <- 0
  }
  dev <- structure(dev, df = choose(p, 2)-m, nobs = ss$n)
  h <- as.numeric(K %*% ss$mean); names(h) <- labels
  tassembly <- proc.time()[[3L]] - t
  structure(list(call = call, graph = graph, decomp = decomp, suffstats = ss,
    precision = K, covariance = Sigma, mean = ss$mean, h = h,
    local_fits = locals, separator_fits = seps, logdet = ld,
    likelihood = list(logLik = ll, trace = trace, deviance = dev),
    convergence = list(converged = TRUE, status = vapply(locals, `[[`, "", "status")),
    diagnostics = list(deviance = dev_message,
      covariance = covariance, control = control),
    timing = list(decomposition = unname(tdec), statistics = unname(tstats),
      local_fitting = unname(tlocal), assembly = unname(tassembly),
      total = unname(proc.time()[[3L]]-start)), engine = engine, version = "0.1.0"),
    class = "atom_gaussian_fit")
}

#' Methods for Gaussian graph fits
#' @param x,object An atom_gaussian_fit.
#' @param ... Unused.
#' @return precision returns the sparse concentration matrix. vcov returns the
#'   fitted variable covariance, or NULL for covariance="none" (not coefficient
#'   uncertainty). logLik returns a logLik object supporting AIC and BIC.
#'   deviance returns the saturated-reference deviance with df and nobs attributes.
#'   summary returns sizes, methods, convergence, likelihood criteria and timing;
#'   print returns its input invisibly.
#' @name atom_gaussian_fit
NULL

#' @rdname atom_gaussian_fit
#' @export
precision <- function(x, ...) UseMethod("precision")
#' @rdname atom_gaussian_fit
#' @export
precision.atom_gaussian_fit <- function(x, ...) x$precision
#' @rdname atom_gaussian_fit
#' @export
vcov.atom_gaussian_fit <- function(object, ...) object$covariance
#' @rdname atom_gaussian_fit
#' @export
logLik.atom_gaussian_fit <- function(object, ...) object$likelihood$logLik
#' @rdname atom_gaussian_fit
#' @export
deviance.atom_gaussian_fit <- function(object, ...) object$likelihood$deviance
#' @rdname atom_gaussian_fit
#' @method print atom_gaussian_fit
#' @export
print.atom_gaussian_fit <- function(x, ...) {
  cat(sprintf("Gaussian graph fit: n=%d, p=%d, edges=%d, atoms=%d\n",
    x$suffstats$n, x$suffstats$p, igraph::ecount(x$graph), length(x$local_fits)))
  cat(sprintf("logLik=%.6f; all local fits successful; total %.3fs\n",
    as.numeric(stats::logLik(x)), x$timing$total))
  invisible(x)
}
#' @rdname atom_gaussian_fit
#' @method summary atom_gaussian_fit
#' @export
summary.atom_gaussian_fit <- function(object, ...) {
  c(list(n = object$suffstats$n), decomp_cost(object$decomp),
    list(incomplete_atom_count = sum(!object$decomp$complete),
      methods = unique(vapply(object$local_fits, function(z) paste(z$engine, z$method, sep = ":"), "")),
      convergence = object$convergence, logLik = stats::logLik(object),
      AIC = stats::AIC(object), BIC = stats::BIC(object), timing = object$timing))
}
