.gaussian_control <- function(control) {
  defaults <- list(tol = 1e-9, contract_tol = 1e-6, maxit = 10000L,
                   method = "ncd", fallback = TRUE)
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      anyNA(names(control)) || any(!nzchar(names(control))) || anyDuplicated(names(control)))))
    stop("control must be a named list without duplicate names.", call. = FALSE)
  if (length(setdiff(names(control), names(defaults)))) stop("Unknown control option.", call. = FALSE)
  defaults[names(control)] <- control
  for (nm in c("tol", "contract_tol", "maxit")) {
    x <- defaults[[nm]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0)
      stop("Control tolerances and maxit must be positive finite scalars.", call. = FALSE)
  }
  if (defaults$maxit != trunc(defaults$maxit) || defaults$maxit > .Machine$integer.max)
    stop("maxit must be a positive integer.", call. = FALSE)
  defaults$method <- match.arg(defaults$method, c("ncd", "covips", "conips"))
  if (!is.logical(defaults$fallback) || length(defaults$fallback) != 1L || is.na(defaults$fallback))
    stop("fallback must be TRUE or FALSE.", call. = FALSE)
  defaults
}

.fit_gaussian_grim <- function(S, graph, n, control) {
  if (!requireNamespace("gRim", quietly = TRUE))
    .gaussian_failure("solver_failed", "Install gRim to fit incomplete atoms with engine='gRim'.")
  e <- t(igraph::as_edgelist(graph, names = FALSE)); storage.mode(e) <- "integer"
  run <- function(method) gRim::fit_ggm_grips(S = S, formula = e, nobs = n,
    method = method, maxit = as.integer(control$maxit), eps = control$tol, convcrit = 1)
  method <- control$method
  z <- run(method)
  fallback <- !isTRUE(z$converged) && method == "ncd" && control$fallback
  if (fallback) { method <- "covips"; z <- run(method) }
  list(K = z$K, Sigma = z$Sigma, iterations = z$details$iter,
       converged = isTRUE(z$converged), method = method,
       message = if (fallback) "NCD did not converge; retried with covips." else "")
}

.fit_gaussian_ggm <- function(S, graph, n, control) {
  if (!requireNamespace("ggm", quietly = TRUE))
    .gaussian_failure("solver_failed", "Install ggm to use engine='ggm'.")
  # fitConGraph starts from S and has no iteration-limit argument.
  .gaussian_chol(S)
  a <- as.matrix(igraph::as_adjacency_matrix(graph, sparse = TRUE))
  dimnames(S) <- dimnames(a) <- list(paste0("V", seq_len(nrow(S))), paste0("V", seq_len(nrow(S))))
  z <- ggm::fitConGraph(amat = a, S = S, n = n, tol = control$tol)
  list(K = chol2inv(.gaussian_chol(z$Shat, "solver_failed")), Sigma = z$Shat,
       iterations = z$it, converged = TRUE, method = "coordinate_descent", message = "")
}

.fit_gaussian_atom <- function(S, graph, n, complete, engine, control,
                               atom_id = 1L, vertices = seq_len(nrow(S))) {
  start <- proc.time()[[3L]]
  result <- tryCatch({
    if (any(diag(S) <= 0)) .gaussian_failure("singular_statistics", "An atom has zero empirical variance.")
    scale <- sqrt(diag(S)); scaling <- outer(scale, scale)
    C <- S/scaling
    if (complete) {
      K <- chol2inv(.gaussian_chol(C))/scaling
      z <- list(K = K, Sigma = S, iterations = 0L, converged = TRUE,
                method = "cholesky", message = "")
      used <- "closed_form"
    } else {
      used <- engine
      # Both APIs fit the supplied covariance without applying n/(n-1).
      # Correlation scaling improves tolerance interpretation and is undone here.
      z <- if (used == "gRim") .fit_gaussian_grim(C, graph, n, control) else .fit_gaussian_ggm(C, graph, n, control)
      if (!isTRUE(z$converged)) .gaussian_failure("nonconverged", "Local solver did not converge; MLE existence is undetermined.")
      q <- nrow(S)
      if (!is.matrix(z$K) || !identical(dim(z$K), c(q, q)) || any(!is.finite(z$K)))
        .gaussian_failure("solver_failed", "Engine returned an invalid precision matrix.")
      if (!is.matrix(z$Sigma) || !identical(dim(z$Sigma), c(q, q)) || any(!is.finite(z$Sigma)))
        .gaussian_failure("solver_failed", "Engine returned an invalid covariance matrix.")
      adj <- as.matrix(igraph::as_adjacency_matrix(graph, sparse = TRUE)) != 0
      diag(adj) <- TRUE
      if (max(abs(z$K-t(z$K))) > control$contract_tol * max(1, abs(z$K)))
        .gaussian_failure("nonconverged", "Local precision violates symmetry or graph zeros.")
      if (max(abs(z$Sigma-t(z$Sigma))) > control$contract_tol * max(1, abs(z$Sigma)))
        .gaussian_failure("solver_failed", "Local covariance violates symmetry.")
      .gaussian_chol(z$Sigma, "solver_failed")
      z$K <- (z$K + t(z$K))/2
      .gaussian_chol(z$K, "solver_failed")
      nonedge_scale <- outer(sqrt(diag(z$K)), sqrt(diag(z$K)))
      if (max(c(0, abs(z$K[!adj])/nonedge_scale[!adj])) > control$contract_tol)
        .gaussian_failure("nonconverged", "Local precision violates symmetry or graph zeros.")
      z$K[!adj] <- 0
      .gaussian_chol(z$K, "solver_failed")
      if (max(abs(z$K %*% z$Sigma - diag(q))) > control$contract_tol)
        .gaussian_failure("solver_failed", "Local precision and covariance matrices are numerically inconsistent.")
      if (max(abs((z$Sigma-C)[adj])) > control$contract_tol)
        .gaussian_failure("nonconverged", "Local fit violates the covariance score equations.")
      z$K <- z$K/scaling; z$Sigma <- z$Sigma*scaling
    }
    dimnames(z$K) <- dimnames(z$Sigma) <- dimnames(S)
    c(list(atom_id = atom_id, vertices = vertices, complete = complete, engine = used),
      z, list(status = if (complete) "closed_form" else "converged",
              elapsed = unname(proc.time()[[3L]]-start)))
  }, error = function(e) {
    if (!inherits(e, "atom_gaussian_error"))
      .gaussian_failure("solver_failed", conditionMessage(e), atom_id, vertices)
    e$atom_id <- atom_id; e$vertices <- vertices
    stop(e)
  })
  result
}

.fit_gaussian_separator <- function(S, vertices, separator_id) {
  z <- tryCatch(chol2inv(.gaussian_chol(S)), error = function(e) {
    e$separator_id <- separator_id; e$vertices <- vertices; stop(e)
  })
  dimnames(z) <- dimnames(S)
  list(separator_id = separator_id, vertices = vertices, Sigma = S, K = z,
       status = "closed_form")
}
