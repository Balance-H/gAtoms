# Audit-only tests for the G-Wishart module: parameterization and the exact
# complete-graph log normalizing constant.
#
# Independence policy: every reference quantity below is reconstructed from
# scratch (explicit arithmetic, a recurrence-based multivariate gamma, a
# Cholesky determinant computed inline, or a direct BDgraph::gnorm() call with
# an independently constructed strict upper-triangular adjacency matrix).  The
# module's own helpers are never used as their own reference.

gwaudit_ns <- function() asNamespace("gAtoms")

gwaudit_graph <- function(kind) {
  g <- switch(kind,
    singleton = igraph::make_empty_graph(1, directed = FALSE),
    edge = igraph::make_full_graph(2),
    complete3 = igraph::make_full_graph(3),
    complete5 = igraph::make_full_graph(5),
    path4 = igraph::make_ring(4, circular = FALSE))
  igraph::set_vertex_attr(g, "name", value = LETTERS[seq_len(igraph::vcount(g))])
}

# Non-diagonal SPD matrix; deliberately not diagonal so that block/order bugs
# cannot hide.
gwaudit_spd <- function(q, seed = 1L, off = 0.25) {
  set.seed(seed)
  A <- matrix(off, q, q)
  diag(A) <- 1
  A + diag(q) * 0.7
}

# Independent exact log multivariate gamma.  Whenever 2a is an integer >= q,
# every log-gamma argument a + (1-j)/2 is a half-integer and lgamma() is exact,
# so this closed form is a rigorous reference:
#   Gamma_q(a) = pi^{q(q-1)/4} prod_{j=1..q} Gamma(a + (1-j)/2)
gwaudit_log_mgamma_exact <- function(a, q) {
  (q * (q - 1) / 4) * log(pi) + sum(lgamma(a + (1 - seq_len(q)) / 2))
}

# Independent check of the defining recurrence
#   Gamma_q(a + 1) / Gamma_q(a) = prod_{j=1..q} (a + (1-j)/2)
# which holds for every a and q.
gwaudit_mgamma_ratio <- function(a, q) prod(a + (1 - seq_len(q)) / 2)

# Independent complete-graph log normalizing constant written out explicitly.
gwaudit_log_norm_complete <- function(D, b) {
  q <- nrow(D)
  a <- (b + q - 1) / 2
  log_det <- 2 * sum(log(diag(chol(D))))
  q * (b + q - 1) / 2 * log(2) + gwaudit_log_mgamma_exact(a, q) - a * log_det
}

# Independently constructed strict upper-triangular adjacency of an induced graph.
gwaudit_upper_adj <- function(g, idx) {
  sub <- igraph::induced_subgraph(g, idx)
  el <- igraph::as_edgelist(sub, names = FALSE)
  q <- igraph::vcount(sub)
  a <- matrix(0, q, q)
  if (nrow(el)) for (r in seq_len(nrow(el))) {
    i <- min(el[r, ]); j <- max(el[r, ])
    a[i, j] <- 1
  }
  a
}

test_that("The independent multivariate-gamma reference matches known values", {
  # q = 1 reduces to lgamma
  expect_equal(gwaudit_log_mgamma_exact(3.5, 1L), lgamma(3.5), tolerance = 1e-13)
  # q = 2 closed form: Gamma_2(a) = sqrt(pi) Gamma(a) Gamma(a - 1/2)
  expect_equal(gwaudit_log_mgamma_exact(4, 2L),
               log(sqrt(pi)) + lgamma(4) + lgamma(3.5), tolerance = 1e-13)
  # q = 3: pi^{3/2} Gamma(a) Gamma(a-1/2) Gamma(a-1)
  expect_equal(gwaudit_log_mgamma_exact(4, 3L),
               1.5 * log(pi) + lgamma(4) + lgamma(3.5) + lgamma(3), tolerance = 1e-13)
})

test_that("The module's multivariate log-gamma matches the independent reference", {
  gm <- get(".gwishart_log_mgamma", envir = gwaudit_ns())
  # rigorous half-integer reference (2a integer >= q), including non-integer a
  for (q in c(1L, 2L, 3L, 4L, 5L)) for (a in c(q / 2, q / 2 + 1, 6.5, 10)) {
    expect_equal(gm(a, q), gwaudit_log_mgamma_exact(a, q), tolerance = 1e-11)
  }
  # the defining recurrence, valid for arbitrary real a
  for (q in c(2L, 3L, 5L)) for (a in c(2.5, 3.7, 4, 9)) {
    expect_equal(exp(gm(a + 1, q) - gm(a, q)), gwaudit_mgamma_ratio(a, q),
                 tolerance = 1e-11)
  }
  # q = 1 is exactly lgamma for arbitrary a
  for (a in c(1.5, 3.7, 12)) expect_equal(gm(a, 1L), lgamma(a), tolerance = 1e-13)
  # the boundary a = (q-1)/2 is a genuine pole of the multivariate gamma
  expect_true(is.infinite(gm(2, 5L)))
  expect_true(is.finite(gm(2.5, 5L)))
})

test_that("The module's Cholesky log-determinant is stable and correct", {
  ldet <- get(".gwishart_log_det", envir = gwaudit_ns())
  D <- gwaudit_spd(4, seed = 5)
  expect_equal(ldet(D), 2 * sum(log(diag(chol(D)))), tolerance = 1e-15)
  expect_equal(ldet(D), log(det(D)), tolerance = 1e-12)
  # scaled matrices must not overflow or lose the log scale
  expect_true(is.finite(ldet(D * 1e6)))
  expect_true(is.finite(ldet(D * 1e-6)))
  expect_equal(ldet(D * 1e6) - ldet(D), 4 * log(1e6), tolerance = 1e-9)
})

test_that("Complete-graph closed form equals the independent explicit formula", {
  for (q in c(1L, 2L, 3L, 5L)) for (b in c(3, 4, 7)) {
    D <- gwaudit_spd(q, seed = q)
    fit <- gwishart_log_norm(gwaudit_graph(switch(as.character(q),
      "1" = "singleton", "2" = "edge", "3" = "complete3", "5" = "complete5")),
      b = b, D = D)
    expect_equal(fit$log_norm, gwaudit_log_norm_complete(D, b), tolerance = 1e-11)
    expect_false(fit$approximate)
    expect_identical(fit$diagnostics$n_bdgraph_calls, 0L)
    expect_true(all(vapply(fit$local_fits, function(z) z$engine == "closed_form", logical(1))))
    expect_true(all(vapply(fit$local_fits, function(z) !isTRUE(z$approximate), logical(1))))
  }
})

test_that("Singleton reduces to the scalar gamma integral", {
  D <- matrix(3.25, 1, 1)
  fit <- gwishart_log_norm(gwaudit_graph("singleton"), b = 4, D = D)
  expect_equal(fit$log_norm, log(2^(4 / 2)) + lgamma(4 / 2) - (4 / 2) * log(3.25),
               tolerance = 1e-13)
  # independent numerical integration of the scalar integral
  # I = int_0^inf k^{(b-2)/2} exp(-k d/2) dk = Gamma(b/2) (d/2)^{-b/2}
  d <- 3.25; bb <- 4
  num <- stats::integrate(function(k) k^((bb - 2) / 2) * exp(-k * d / 2),
                          0, Inf, rel.tol = 1e-13)$value
  expect_equal(fit$log_norm, log(num), tolerance = 1e-10)
})

test_that("Complete-graph values agree with the complete branch of BDgraph::gnorm", {
  skip_if_not_installed("BDgraph")
  for (q in c(2L, 3L, 5L)) for (b in c(3, 4, 6)) {
    g <- gwaudit_graph(switch(as.character(q), "2" = "edge", "3" = "complete3",
                              "5" = "complete5"))
    D <- gwaudit_spd(q, seed = 100 + q)
    fit <- gwishart_log_norm(g, b = b, D = D)
    ref <- BDgraph::gnorm(gwaudit_upper_adj(g, seq_len(q)), b = b, D = D, iter = 100L)
    expect_equal(fit$log_norm, ref, tolerance = 1e-10)
  }
})

test_that("The complete-graph constant depends on D only through log|D| (theory check)", {
  q <- 3L; b <- 4
  D <- gwaudit_spd(q, seed = 7)
  g <- gwaudit_graph("complete3")
  a <- (b + q - 1) / 2
  # scaling D by c must shift the constant by -a q log c (log|D| -> log|D| + q log c)
  for (cc in c(0.5, 2, 10)) {
    expect_equal(gwishart_log_norm(g, b = b, D = cc * D)$log_norm,
                 gwishart_log_norm(g, b = b, D = D)$log_norm - a * q * log(cc),
                 tolerance = 1e-10)
  }
})

test_that("A structurally non-diagonal but determinant-equivalent D gives the same value", {
  # for a complete graph only log|D| and q enter, so any SPD D with the same
  # determinant must produce the same constant; this is an independent check
  # that no off-diagonal or per-entry structure leaks into the formula
  D1 <- matrix(c(2, 0.6, 0.6, 2), 2, 2)
  D2 <- diag(c(4, det(D1) / 4))
  expect_equal(det(D1), det(D2), tolerance = 1e-14)
  g <- gwaudit_graph("edge")
  expect_equal(gwishart_log_norm(g, b = 3, D = D1)$log_norm,
               gwishart_log_norm(g, b = 3, D = D2)$log_norm, tolerance = 1e-11)
})

test_that("Complete-prime and separator terms use the same closed-form convention", {
  # separators are complete cliques: their stored value must equal the explicit
  # complete-graph formula applied to the corresponding principal D block
  g <- gwaudit_graph("path4")
  D <- gwaudit_spd(4, seed = 9)
  fit <- gwishart_log_norm(g, b = 5, D = D)
  expect_gt(length(fit$separator_fits), 0L)
  for (z in fit$separator_fits) {
    idx <- match(z$vertices, igraph::V(fit$graph)$name)
    expect_identical(z$engine, "closed_form")
    expect_false(isTRUE(z$approximate))
    expect_equal(z$log_norm, gwaudit_log_norm_complete(D[idx, idx, drop = FALSE], 5),
                 tolerance = 1e-12)
  }
  for (z in fit$local_fits) {
    idx <- match(z$vertices, igraph::V(fit$graph)$name)
    expect_equal(z$log_norm, gwaudit_log_norm_complete(D[idx, idx, drop = FALSE], 5),
                 tolerance = 1e-12)
  }
})
