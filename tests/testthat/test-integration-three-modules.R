# Package-level integration tests: public API visibility and reuse of one
# atom_decomp object by all three statistical modules.
#
# These tests are architectural: they do not re-check the statistics already
# covered by the Gaussian, log-linear and G-Wishart audit test files. They
# verify that the integrated package exposes the intended API and that a single
# decomposition object can be shared by the Gaussian, log-linear and G-Wishart
# modules without mutation.

integration_graph <- function() {
  igraph::set_vertex_attr(igraph::make_ring(4, circular = FALSE),
                          "name", value = c("A", "B", "C", "D"))
}

integration_cycle <- function() {
  igraph::set_vertex_attr(igraph::make_ring(4),
                          "name", value = c("A", "B", "C", "D"))
}

test_that("The public API of the three statistical modules is visible", {
  for (fn in c("atom_decompose", "atoms", "separators", "residuals",
               "prime_ordering", "atom_tree", "validate", "decomp_cost",
               "decompose_atoms", "get_minimal_collapsible",
               "as_igraph", "as_dgCMatrix",
               "gaussian_suffstats", "fit_gaussian_graph", "precision",
               "margin_counts", "fit_loglinear_graph",
               "gwishart_log_norm")) {
    expect_true(exists(fn, mode = "function"), info = fn)
  }
  expect_true(is.function(gwishart_log_norm))
  expect_identical(names(formals(gwishart_log_norm)),
                   c("graph", "b", "D", "decomp", "engine", "control"))
  # the legacy Gaussian interface and the removed plot convenience wrapper
  # must remain absent, while the plot method stays available
  expect_false(exists("fit_ggm_localized", mode = "function"))
  expect_false(exists("plot_atom_tree", mode = "function"))
  expect_true("plot.atom_decomp" %in% as.character(utils::methods("plot")))
})

test_that("One atom_decomp is reused by all three modules without mutation", {
  G <- integration_graph()
  dec <- atom_decompose(G)
  dec_before <- list(
    atoms = atoms(dec),
    separators = separators(dec),
    residuals = residuals(dec),
    prime_ordering = prime_ordering(dec),
    vertex_map = dec$vertex_map,
    components = dec$components,
    parent = dec$parent,
    children = dec$children,
    complete = dec$complete
  )

  # --- Gaussian: a path graph has only complete prime components, so the
  #     closed-form path is used and no optional engine is required.
  set.seed(11)
  X <- matrix(stats::rnorm(40 * 4), nrow = 40, ncol = 4)
  colnames(X) <- c("A", "B", "C", "D")
  fit_g <- fit_gaussian_graph(X, G, decomp = dec)

  # --- log-linear: all prime components are complete (closed form only)
  dat <- expand.grid(A = c("a", "b"), B = c("a", "b"),
                     C = c("a", "b"), D = c("a", "b"))
  dat <- dat[rep(seq_len(nrow(dat)), 1L + (seq_len(nrow(dat)) * 5L) %% 7L), , drop = FALSE]
  fit_l <- fit_loglinear_graph(dat, G, decomp = dec, engine = "loglin")

  # --- G-Wishart: complete primes only, hence exact and BDgraph-free
  D <- diag(4)
  fit_w <- gwishart_log_norm(G, b = 3, D = D, decomp = dec)

  # expected classes
  expect_s3_class(fit_g, "atom_gaussian_fit")
  expect_s3_class(fit_l, "atom_loglinear_fit")
  expect_s3_class(fit_w, "atom_gwishart_norm")

  # finite primary quantities
  expect_true(is.finite(as.numeric(stats::logLik(fit_g))))
  expect_true(is.finite(as.numeric(stats::logLik(fit_l))))
  expect_true(is.finite(fit_w$log_norm))

  # Gaussian precision is a finite 4x4 sparse symmetric matrix
  prec <- precision(fit_g)
  expect_true(all(dim(as.matrix(prec)) == c(4L, 4L)))
  expect_true(all(is.finite(as.matrix(prec))))

  # log-linear deviance and factorized fitted output
  expect_true(is.finite(stats::deviance(fit_l)))
  fac <- stats::fitted(fit_l)
  expect_true(is.list(fac))
  expect_length(fac$atoms, length(fit_l$local_fits))
  expect_identical(names(fac$levels), c("A", "B", "C", "D"))

  # G-Wishart exactness bookkeeping
  expect_false(fit_w$approximate)
  expect_identical(fit_w$engine, "closed_form")
  expect_identical(fit_w$diagnostics$n_bdgraph_calls, 0L)
  expect_true(all(is.finite(fit_w$atom_log_norm)))
  expect_true(all(is.finite(fit_w$separator_log_norm)))
  expect_equal(fit_w$log_norm,
               sum(fit_w$atom_log_norm) - sum(fit_w$separator_log_norm),
               tolerance = 1e-15)

  # decomposition reuse and non-mutation
  expect_true(fit_w$diagnostics$decomp_reused)
  expect_identical(atoms(dec), dec_before$atoms)
  expect_identical(separators(dec), dec_before$separators)
  expect_identical(residuals(dec), dec_before$residuals)
  expect_identical(prime_ordering(dec), dec_before$prime_ordering)
  expect_identical(dec$vertex_map, dec_before$vertex_map)
  expect_identical(dec$components, dec_before$components)
  expect_identical(dec$parent, dec_before$parent)
  expect_identical(dec$children, dec_before$children)
  expect_identical(dec$complete, dec_before$complete)
  expect_true(isTRUE(validate(dec, expensive = TRUE)))

  # every module keeps the same graph vertex labels
  expect_identical(igraph::V(fit_g$graph)$name, c("A", "B", "C", "D"))
  expect_identical(igraph::V(fit_l$graph)$name, c("A", "B", "C", "D"))
  expect_identical(igraph::V(fit_w$graph)$name, c("A", "B", "C", "D"))
  expect_identical(dec$vertex_map$label, c("A", "B", "C", "D"))
})

test_that("Supplying a decomposition avoids a second decomposition", {
  G <- integration_graph()
  dec <- atom_decompose(G)
  fit_w <- gwishart_log_norm(G, b = 3, D = diag(4), decomp = dec)
  expect_identical(fit_w$timing$decomposition, 0)
  # a decomposition of a different graph is rejected
  expect_error(gwishart_log_norm(integration_cycle(), b = 3, D = diag(4),
                                 decomp = dec))
})

test_that("G-Wishart closed-form case works without BDgraph", {
  G <- integration_graph()
  fit <- gwishart_log_norm(G, b = 3, D = diag(4))
  expect_identical(fit$diagnostics$n_bdgraph_calls, 0L)
  expect_false(fit$approximate)
  expect_true(is.finite(fit$log_norm))
})

test_that("G-Wishart incomplete-prime case uses BDgraph when available", {
  skip_if_not_installed("BDgraph")
  G <- integration_cycle()
  set.seed(1)
  fit <- gwishart_log_norm(G, b = 3, D = diag(4), control = list(iter = 50L))
  expect_true(fit$approximate)
  expect_identical(fit$engine, "decomposition+BDgraph")
  expect_identical(fit$diagnostics$n_bdgraph_calls, 1L)
  expect_true(is.finite(fit$log_norm))
})

test_that("Optional Gaussian and log-linear engines stay optional", {
  G <- integration_graph()
  set.seed(3)
  X <- matrix(stats::rnorm(30 * 4), nrow = 30, ncol = 4)
  colnames(X) <- c("A", "B", "C", "D")
  expect_s3_class(fit_gaussian_graph(X, G), "atom_gaussian_fit")
  if (requireNamespace("gRim", quietly = TRUE))
    expect_s3_class(fit_gaussian_graph(X, G, engine = "gRim"), "atom_gaussian_fit")
  if (requireNamespace("ggm", quietly = TRUE))
    expect_s3_class(fit_gaussian_graph(X, G, engine = "ggm"), "atom_gaussian_fit")

  dat <- expand.grid(A = c("a", "b"), B = c("a", "b"),
                     C = c("a", "b"), D = c("a", "b"))
  dat <- dat[rep(seq_len(nrow(dat)), 1L + (seq_len(nrow(dat)) * 3L) %% 5L), , drop = FALSE]
  expect_s3_class(fit_loglinear_graph(dat, G, engine = "loglin"),
                  "atom_loglinear_fit")
  if (requireNamespace("gRim", quietly = TRUE))
    expect_s3_class(fit_loglinear_graph(dat, G, engine = "gRim"),
                    "atom_loglinear_fit")
})
