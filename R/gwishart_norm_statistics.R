# gAtoms G-Wishart normalizing constant -- statistical utilities.
#
# Parameterization (identical to BDgraph::gnorm):
#
#   p(K | G, b, D) = I_G(b, D)^-1 |K|^{(b-2)/2} exp(-tr(K D)/2),   K in P_G
#
# with P_G the cone of symmetric positive definite matrices with zeros where G
# has no edge.  This file provides only input validation and the exact
# complete-graph log normalizing constant.  No Monte Carlo work happens here.

# ---------------------------------------------------------------------------
# Validation of the G-Wishart parameters b and D
# ---------------------------------------------------------------------------

# Shape parameter b, restricted to integer values b >= 3.
#
# The density needs b > 2 for the exponent (b-2)/2 to be usable on the whole
# cone, and the complete-graph closed form below is in fact valid for every
# b > 0 (all log-gamma arguments stay positive).  BDgraph::gnorm(), however,
# uses an integer-valued b in its Monte Carlo backend for incomplete graphs:
# that branch passes b to its compiled kernel through integer coercion while
# its R-level bookkeeping term uses the real value.  Requiring integer b
# globally keeps complete-prime, separator and incomplete-prime contributions
# under one coherent G-Wishart parameterization, avoids depending on that
# mixed interface, and makes the public contract independent of whether the
# current graph happens to contain an incomplete prime.  Values that cannot be
# represented as an R integer are rejected as well so that no downstream
# overflow can occur.  Invalid input is never repaired by rounding.
.gwishart_validate_b <- function(b) {
  if (!is.numeric(b) || length(b) != 1L || !is.finite(b) || b < 3 || b != round(b))
    stop("b must be a finite integer-valued scalar with b >= 3.", call. = FALSE)
  if (b > .Machine$integer.max)
    stop("b exceeds the integer range supported by the BDgraph backend.",
         call. = FALSE)
  as.numeric(b)
}

# Scale matrix D.  Exact requirements only: numeric, square, dimensionally
# compatible, symmetric, finite and positive definite.  Nothing is regularized
# and nothing is repaired.
.gwishart_validate_D <- function(D, q) {
  if (is.data.frame(D)) D <- as.matrix(D)
  if (!is.matrix(D) || !is.numeric(D))
    stop("D must be a numeric matrix.", call. = FALSE)
  if (nrow(D) != ncol(D))
    stop("D must be square.", call. = FALSE)
  if (nrow(D) != q)
    stop(sprintf("D has dimension %d but the graph has %d vertices.", nrow(D), q),
         call. = FALSE)
  if (anyNA(D) || any(!is.finite(D)))
    stop("D must not contain NA, NaN or infinite values.", call. = FALSE)
  if (!isTRUE(all.equal(D, t(D), tolerance = 1e-8, check.attributes = FALSE)))
    stop("D must be symmetric.", call. = FALSE)
  D <- (D + t(D)) / 2
  if (any(eigen(D, symmetric = TRUE, only.values = TRUE)$values <= 0))
    stop("D must be positive definite.", call. = FALSE)
  D
}

# ---------------------------------------------------------------------------
# Stable log determinant
# ---------------------------------------------------------------------------

# log|A| through the Cholesky factor of a symmetric positive definite A, i.e.
# 2 * sum(log(diag(chol(A)))).  det()/log(det()) is deliberately never used.
.gwishart_log_det <- function(A) {
  R <- chol(A)
  2 * sum(log(diag(R)))
}

# ---------------------------------------------------------------------------
# Multivariate log-gamma and the exact complete-graph constant
# ---------------------------------------------------------------------------

# log Gamma_q(a) = q(q-1)/4 log(pi) + sum_{j=1..q} log Gamma(a + (1-j)/2)
.gwishart_log_mgamma <- function(a, q) {
  if (q < 1L) stop("q must be a positive integer.", call. = FALSE)
  (q * (q - 1) / 4) * log(pi) + sum(lgamma(a + (1 - seq_len(q)) / 2))
}

# log I_{K_q}(b, D) for a complete graph on q vertices:
#
#   I_{K_q}(b, D) = 2^{q(b+q-1)/2} Gamma_q((b+q-1)/2) |D|^{-(b+q-1)/2}
#
# hence, entirely in log scale,
#
#   log I_{K_q}(b, D) = q(b+q-1)/2 log 2
#                     + log Gamma_q((b+q-1)/2)
#                     - (b+q-1)/2 log|D|.
#
# D must already be ordered according to the local vertex order.  The q = 1
# case reduces to the scalar gamma integral
# log I = log 2^{b/2} + lgamma(b/2) - (b/2) log D_{11}.
.gwishart_log_norm_complete <- function(D, b) {
  if (!is.matrix(D) || nrow(D) != ncol(D))
    stop("D must be a square numeric matrix.", call. = FALSE)
  q <- nrow(D)
  if (q < 1L) stop("D must have at least one row and column.", call. = FALSE)
  b <- .gwishart_validate_b(b)
  D <- .gwishart_validate_D(D, q)
  a <- (b + q - 1) / 2
  (q * (b + q - 1) / 2) * log(2) + .gwishart_log_mgamma(a, q) - a * .gwishart_log_det(D)
}

# ---------------------------------------------------------------------------
# Vertex ordering of D
# ---------------------------------------------------------------------------

# Rows/columns of D must correspond to the canonical vertex order of the
# graph.  A named D is reordered deterministically by name; an unnamed D is
# taken to be already in canonical order and is never reordered.
.gwishart_align_D <- function(D, variables, p) {
  if (is.data.frame(D)) D <- as.matrix(D)
  if (!is.matrix(D) || !is.numeric(D))
    stop("D must be a numeric matrix.", call. = FALSE)
  if (nrow(D) != ncol(D))
    stop("D must be square.", call. = FALSE)
  rn <- rownames(D); cn <- colnames(D)
  named <- !is.null(rn) || !is.null(cn)
  if (named) {
    if (is.null(rn) || is.null(cn))
      stop("D must have both row and column names, or neither.", call. = FALSE)
    if (!identical(rn, cn))
      stop("D row names and column names must agree.", call. = FALSE)
    if (anyNA(rn) || anyDuplicated(rn))
      stop("D row and column names must be unique and nonmissing.", call. = FALSE)
    if (!setequal(rn, variables))
      stop("D row names must match the graph vertex names exactly.", call. = FALSE)
    D <- D[variables, variables, drop = FALSE]
  }
  D <- .gwishart_validate_D(D, p)
  D
}

# ---------------------------------------------------------------------------
# Local quantities on a vertex subset
# ---------------------------------------------------------------------------

# Principal submatrix D[A, A] in the local vertex order.
.gwishart_block <- function(D, index) {
  D[index, index, drop = FALSE]
}

# Binary symmetric adjacency matrix of the induced subgraph on `index`, in the
# same local order as the corresponding principal block of D.  Note that
# igraph::as_adjacency_matrix() only fills the upper triangle for an undirected
# graph, so the result is explicitly symmetrized here; the strict
# upper-triangular form required by BDgraph is produced separately in
# gwishart_norm_engine.R.
.gwishart_adjacency <- function(graph, index) {
  sub <- igraph::induced_subgraph(graph, index)
  a <- igraph::as_adjacency_matrix(sub, sparse = FALSE, type = "both")
  q <- nrow(a)
  a <- matrix(as.numeric(a != 0), nrow = q, ncol = q,
              dimnames = list(igraph::V(sub)$name, igraph::V(sub)$name))
  a[lower.tri(a)] <- t(a)[lower.tri(a)]
  diag(a) <- 0
  a
}
