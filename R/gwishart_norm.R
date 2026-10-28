# gAtoms G-Wishart normalizing constant -- decomposition-aware orchestration.
#
# Mathematical mapping
# --------------------
# For a perfect sequence of prime components P_1, ..., P_k with clique
# separators S_2, ..., S_k (the ordered separators produced by
# gAtoms::atom_decompose()), the G-Wishart normalizing constant factorizes as
#
#   I_G(b, D) = prod_i I_{G_{P_i}}(b, D_{P_i}) / prod_{i>=2} I_{G_{S_i}}(b, D_{S_i})
#
# and therefore, entirely in log scale,
#
#   log I_G(b, D) = sum_i log I_{G_{P_i}}(b, D_{P_i})
#                 - sum_{i>=2} log I_{G_{S_i}}(b, D_{S_i}),
#
# with D_A = D[A, A] and A a vertex subset in the local vertex order.
#
# Responsibilities kept strictly separated:
#   * the decomposition belongs to gAtoms (atom_decompose()) and is reused,
#     never recomputed when the caller supplies it and never reimplemented;
#   * separators are cliques by construction and are always evaluated with the
#     exact complete-graph closed form;
#   * complete prime components are evaluated with the same exact closed form;
#   * incomplete prime components are the ONLY terms handed to BDgraph::gnorm().
#
# The decomposition identity itself is exact.  When every prime component is
# complete the whole result is exact; when incomplete primes are present the
# global value is
#
#   sum_i log I_{P_i} - sum_j log I_{S_j}
#
# with exact separator and complete-prime terms and Monte Carlo approximations
# only for the incomplete-prime terms.  The result object therefore never
# claims exactness once BDgraph has been used (approximate = TRUE).
#
# Disconnected graphs are handled componentwise: I_G = prod_c I_{G_c} and no
# artificial separators are created between connected components.  Only the
# principal block D[A, A] of a local vertex set A enters a computation, so
# off-block entries of D are never used.

# Named list of recognized control options and their defaults.
.gwishart_control <- function(control) {
  defaults <- list(iter = 100L)
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      anyNA(names(control)) || anyDuplicated(names(control)) ||
      any(!names(control) %in% names(defaults)))))
    stop("control must be a named list of recognized options (iter).", call. = FALSE)
  defaults[names(control)] <- control
  iter <- defaults$iter
  if (!is.numeric(iter) || length(iter) != 1L || !is.finite(iter) ||
      iter < 1 || iter != trunc(iter) || iter > .Machine$integer.max)
    stop("control$iter must be a single finite integer >= 1.", call. = FALSE)
  list(iter = as.integer(iter))
}

# Reuse a caller-supplied decomposition, or compute it exactly once.  The
# supplied object is checked for basic graph/vertex compatibility and through
# the existing gAtoms validator; no second decomposition algorithm and no
# second structural validator are introduced here.
#
# Timing contract: elapsed is the wall-clock time of the atom_decompose() call
# itself when the decomposition has to be computed, and exactly 0 when the
# caller supplied one.  No other step performed here (class check, validator,
# graph/vertex-order comparison) contributes to it; those costs still fall
# inside the caller's overall total.
.gwishart_decomposition <- function(graph, decomp, variables, p) {
  if (is.null(decomp)) {
    t0 <- proc.time()[[3L]]
    decomp <- atom_decompose(graph)
    return(list(decomp = decomp, reused = FALSE,
                elapsed = unname(proc.time()[[3L]] - t0)))
  }
  if (!inherits(decomp, "atom_decomp"))
    stop("decomp must be an atom_decomp object.", call. = FALSE)
  valid <- tryCatch(validate(decomp, expensive = TRUE), error = function(e) e)
  if (inherits(valid, "error"))
    stop("Invalid decomposition: ", conditionMessage(valid), call. = FALSE)
  if (!isTRUE(valid))
    stop("Invalid decomposition: ", paste(valid$errors, collapse = "; "), call. = FALSE)
  dcmp <- decomp$components$vertices
  if (is.null(dcmp) || length(dcmp) != p)
    stop("decomp does not cover the graph vertices.", call. = FALSE)
  if (igraph::vcount(decomp$graph) != p || !identical(igraph::V(decomp$graph)$name, variables))
    stop("decomp must match graph vertex names and order.", call. = FALSE)
  if (!isTRUE(all.equal(as.numeric(as.matrix(as_dgCMatrix(graph))),
                        as.numeric(as.matrix(as_dgCMatrix(decomp$graph))),
                        tolerance = 0, check.attributes = FALSE)))
    stop("decomp must match graph edges.", call. = FALSE)
  list(decomp = decomp, reused = TRUE, elapsed = 0)
}

# Compact per-component summary of the local contributions.
.gwishart_local_summary <- function(local_fits, separator_fits) {
  if (length(local_fits)) {
    atoms <- data.frame(
      kind = "atom",
      id = vapply(local_fits, function(z) z$atom_id, 0L),
      component_id = vapply(local_fits, function(z) z$component_id, 0L),
      size = vapply(local_fits, function(z) z$size, 0L),
      complete = vapply(local_fits, function(z) isTRUE(z$complete), logical(1)),
      engine = vapply(local_fits, function(z) z$engine, ""),
      approximate = vapply(local_fits, function(z) isTRUE(z$approximate), logical(1)),
      log_norm = vapply(local_fits, function(z) z$log_norm, 0.0),
      elapsed = vapply(local_fits, function(z) z$elapsed, 0.0),
      stringsAsFactors = FALSE)
  } else {
    atoms <- data.frame(kind = character(), id = integer(), component_id = integer(),
      size = integer(), complete = logical(), engine = character(),
      approximate = logical(), log_norm = numeric(), elapsed = numeric(),
      stringsAsFactors = FALSE)
  }
  if (length(separator_fits)) {
    separators <- data.frame(
      kind = "separator",
      id = vapply(separator_fits, function(z) z$separator_id, 0L),
      component_id = vapply(separator_fits, function(z) z$component_id, 0L),
      size = vapply(separator_fits, function(z) z$size, 0L),
      complete = TRUE,
      engine = vapply(separator_fits, function(z) z$engine, ""),
      approximate = FALSE,
      log_norm = vapply(separator_fits, function(z) z$log_norm, 0.0),
      elapsed = vapply(separator_fits, function(z) z$elapsed, 0.0),
      stringsAsFactors = FALSE)
  } else {
    separators <- atoms[0L, , drop = FALSE]
  }
  rbind(atoms, separators)
}

#' Decomposition-aware G-Wishart log normalizing constant
#'
#' Computes \eqn{\log I_G(b, D)} for the G-Wishart distribution on a fixed
#' undirected graph \eqn{G} under the parameterization
#'
#' \deqn{p(K \mid G, b, D) = I_G(b, D)^{-1} |K|^{(b - 2)/2}
#'   \exp\left\{-\frac{1}{2}\operatorname{tr}(KD)\right\},}{%
#'   p(K | G, b, D) = I_G(b, D)^-1 |K|^((b-2)/2) exp(-tr(K D)/2),}
#'
#' for \eqn{K} in the cone of symmetric positive definite matrices with zeros
#' wherever \eqn{G} has no edge. This is the same parameterization used by
#' \code{\link[BDgraph]{gnorm}}, so the value is comparable with it.
#'
#' The computation exploits the atom / prime decomposition of the graph:
#' complete prime components and clique separators are evaluated with the exact
#' complete-graph closed form
#' \eqn{2^{q(b+q-1)/2} \Gamma_q((b+q-1)/2) |D|^{-(b+q-1)/2}}{%
#' 2^(q(b+q-1)/2) Gamma_q((b+q-1)/2) |D|^-(b+q-1)/2}, while incomplete prime
#' components are the only terms passed to \code{BDgraph::gnorm()}. The
#' decomposition identity
#' \eqn{\log I_G = \sum_i \log I_{P_i} - \sum_j \log I_{S_j}}{%
#' log I_G = sum_i log I_Pi - sum_j log I_Sj} is exact, so the returned value is
#' exact when every prime component is complete and approximate only when an
#' incomplete prime requires BDgraph Monte Carlo. Disconnected graphs are
#' handled componentwise with no artificial separators.
#'
#' @param graph An undirected simple \code{igraph}, a binary adjacency matrix or
#'   a \code{dgCMatrix}. Vertex names define the variables of \code{D}. Use the
#'   same graph object for the companion Gaussian and log-linear modules.
#' @param b G-Wishart shape parameter. Must be a finite integer-valued scalar
#'   with \eqn{b \ge 3} and \eqn{b \le} \code{.Machine$integer$max}; the integer
#'   domain is required because the incomplete-prime backend uses an
#'   integer-valued \code{b}, and every local term must share one coherent
#'   parameterization.
#' @param D Numeric symmetric positive definite scale matrix with as many rows
#'   and columns as the graph has vertices. If \code{D} has row and column
#'   names, the names must agree, be unique and match the graph vertices, and
#'   \code{D} is reordered to the canonical graph vertex order; an unnamed
#'   \code{D} is taken to be already in that order and is never reordered.
#' @param decomp Optional compatible \code{atom_decomp} from
#'   \code{\link{atom_decompose}}. It is validated and reused, never recomputed.
#' @param engine Character; the incomplete-prime backend. \code{"BDgraph"} uses
#'   \code{BDgraph::gnorm()} and is the only currently supported value.
#' @param control Named list of control options: \code{iter} (default
#'   \code{100L}) is the number of Monte Carlo iterations passed to
#'   \code{BDgraph::gnorm()}. No other options are recognized.
#'
#' @return An object of class \code{"atom_gwishart_norm"}, a list with elements
#'   \describe{
#'     \item{log_norm}{the value \eqn{\log I_G(b, D)}.}
#'     \item{b, D, graph, decomp}{the validated inputs; \code{D} is stored in
#'       canonical graph vertex order.}
#'     \item{atom_log_norm, separator_log_norm}{named numeric vectors of the
#'       prime and separator contributions.}
#'     \item{local_fits, separator_fits}{per-component records, including the
#'       vertices, size, \code{complete}, \code{engine}, \code{approximate},
#'       \code{log_norm} and the strict upper-triangular \code{adjacency}.}
#'     \item{local_summary}{compact data frame of all local contributions.}
#'     \item{approximate}{logical; \code{TRUE} when at least one incomplete
#'       prime was evaluated by BDgraph Monte Carlo.}
#'     \item{diagnostics}{counts and sizes, including \code{n_atoms},
#'       \code{n_separators}, \code{n_complete_atoms}, \code{n_incomplete_atoms},
#'       \code{max_incomplete_atom_size}, \code{n_bdgraph_calls} and
#'       \code{decomp_reused}.}
#'     \item{timing, control, engine, call}{timing breakdown, the validated
#'       control list, the global engine description
#'       (\code{"closed_form"} or \code{"decomposition+BDgraph"}) and the
#'       matched call.}
#'   }
#'
#' @section Scale matrix handling: \code{D} must be finite, symmetric and
#'   positive definite; it is never repaired and no ridge or jitter is added.
#'   Matrices that are symmetric within a relative tolerance of \code{1e-8} are
#'   internally symmetrized as \code{(D + t(D)) / 2} before the Cholesky
#'   log-determinant is evaluated; clearly asymmetric matrices are rejected.
#'
#' @section Optional dependency: BDgraph is an optional dependency. Graphs whose
#'   prime components are complete need no BDgraph installation; a graph with an
#'   incomplete prime fails with an informative error when BDgraph (or
#'   \code{BDgraph::gnorm}) is not available.
#'
#' @examples
#' G <- igraph::set_vertex_attr(igraph::make_ring(4, circular = FALSE),
#'                             "name", value = c("A", "B", "C", "D"))
#' dec <- atom_decompose(G)
#' D <- diag(4)
#' gwishart_log_norm(G, b = 3, D = D, decomp = dec)
#'
#' @seealso \code{\link{atom_decompose}}, \code{\link{fit_gaussian_graph}},
#'   \code{\link{fit_loglinear_graph}}
#' @export
gwishart_log_norm <- function(graph, b, D, decomp = NULL, engine = "BDgraph",
                              control = list()) {
  start <- proc.time()[[3L]]; call <- match.call()
  engine <- match.arg(engine, c("BDgraph", "closed_form"))
  if (identical(engine, "closed_form"))
    stop("engine = \"closed_form\" cannot evaluate incomplete prime components; ",
         "it has no Monte Carlo backend. Use engine = \"BDgraph\".", call. = FALSE)
  control <- .gwishart_control(control)
  b <- .gwishart_validate_b(b)
  norm <- normalize_graph(graph); graph <- norm$graph
  p <- igraph::vcount(graph)
  if (!p) stop("A G-Wishart normalizing constant requires a nonempty graph.", call. = FALSE)
  variables <- norm$vertex_map$label
  D <- .gwishart_align_D(D, variables, p)
  sel <- .gwishart_decomposition(graph, decomp, variables, p)
  decomp <- sel$decomp
  # timing$decomposition measures only atom_decompose(graph) and is exactly 0
  # when the caller supplied the decomposition; every other step above (engine
  # matching, control/b/D validation, normalize_graph, D alignment and the
  # supplied-decomposition validator) is charged to timing$total instead.
  tdec <- sel$elapsed
  local_start <- proc.time()[[3L]]
  atoms_index <- decomp$atoms_index; separators_index <- decomp$separators_index
  components <- as.integer(decomp$components$vertices)
  k <- length(atoms_index)
  local_fits <- lapply(seq_len(k), function(i) {
    idx <- sort(atoms_index[[i]])
    .gwishart_local_fit(atom_id = i, component_id = components[idx[1L]], index = idx,
      graph = graph, vertices = variables[idx], D = D, b = b,
      complete = isTRUE(decomp$complete[i]), iter = control$iter)
  })
  j <- length(separators_index)
  separator_fits <- lapply(seq_len(j), function(i) {
    idx <- sort(separators_index[[i]])
    .gwishart_separator_fit(separator_id = i, component_id = components[idx[1L]],
      index = idx, vertices = variables[idx], D = D, b = b)
  })
  tlocal <- unname(proc.time()[[3L]] - local_start)
  comp_start <- proc.time()[[3L]]
  if (k) names(local_fits) <- paste0("atom", seq_len(k))
  if (j) names(separator_fits) <- paste0("sep", seq_len(j))
  atom_log_norm <- vapply(local_fits, function(z) z$log_norm, 0.0)
  separator_log_norm <- vapply(separator_fits, function(z) z$log_norm, 0.0)
  if (length(separator_fits))
    names(separator_log_norm) <- vapply(separator_fits, function(z)
      paste0("sep", z$separator_id, ":", paste(z$vertices, collapse = ",")), "")
  log_norm <- sum(atom_log_norm) - sum(separator_log_norm)
  sizes <- lengths(atoms_index)
  incomplete_sizes <- sizes[!as.logical(decomp$complete)]
  n_bdgraph <- sum(vapply(local_fits, function(z) identical(z$engine, "BDgraph"), logical(1)))
  approximate <- n_bdgraph > 0L
  composition <- unname(proc.time()[[3L]] - comp_start)
  structure(list(
    log_norm = log_norm, b = b, D = D, graph = graph, decomp = decomp,
    atom_log_norm = atom_log_norm, separator_log_norm = separator_log_norm,
    local_fits = local_fits, separator_fits = separator_fits,
    local_summary = .gwishart_local_summary(local_fits, separator_fits),
    approximate = approximate,
    diagnostics = list(
      p = p, n_components = length(unique(components)), n_atoms = k,
      n_separators = j, n_complete_atoms = sum(as.logical(decomp$complete)),
      n_incomplete_atoms = k - sum(as.logical(decomp$complete)),
      max_atom_size = if (k) max(sizes) else 0L,
      max_incomplete_atom_size = if (length(incomplete_sizes)) max(incomplete_sizes) else 0L,
      n_bdgraph_calls = n_bdgraph, full_graph_size = p,
      atom_sizes = sizes, separator_sizes = lengths(separators_index),
      decomp_reused = sel$reused,
      note = if (approximate)
        "decomposition exact; incomplete-prime normalizing constants approximated by BDgraph Monte Carlo" else
        "decomposition exact; all prime components and separators evaluated in closed form"),
    timing = list(decomposition = tdec, local = tlocal, composition = composition,
      total = unname(proc.time()[[3L]] - start)),
    control = control,
    engine = if (approximate) "decomposition+BDgraph" else "closed_form",
    call = call), class = "atom_gwishart_norm")
}
