# Audit-only tests for the G-Wishart module: prime/separator factorization,
# disconnected graphs, unit handling and decomposition reuse.

gwaudit_ns <- function() asNamespace("gAtoms")

gwaudit_fun <- function(name) {
  if (exists(name, envir = globalenv(), inherits = FALSE)) return(get(name, envir = globalenv()))
  get(name, envir = gwaudit_ns())
}

gwaudit_graph <- function(kind) {
  g <- switch(kind,
    singleton = igraph::make_empty_graph(1, directed = FALSE),
    edge = igraph::make_full_graph(2),
    path4 = igraph::make_ring(4, circular = FALSE),
    path6 = igraph::make_ring(6, circular = FALSE),
    cliques3 = igraph::graph_from_edgelist(rbind(c(1, 2), c(2, 3), c(1, 3),
                                                 c(3, 4), c(4, 5), c(3, 5)), directed = FALSE),
    cycle4 = igraph::make_ring(4),
    disc_complete = igraph::disjoint_union(igraph::make_full_graph(3),
                                           igraph::make_full_graph(2)),
    disc_mixed = igraph::disjoint_union(igraph::make_full_graph(3),
                                        igraph::make_ring(4)),
    disc_singletons = igraph::make_empty_graph(3, directed = FALSE),
    mixed = igraph::graph_from_edgelist(rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1),
                                              c(3, 5), c(5, 6)), directed = FALSE))
  igraph::set_vertex_attr(g, "name", value = LETTERS[seq_len(igraph::vcount(g))])
}

gwaudit_spd <- function(q, seed = 1L, off = 0.25) {
  set.seed(seed); A <- matrix(off, q, q); diag(A) <- 1; A + diag(q) * 0.7
}

gwaudit_log_mgamma_exact <- function(a, q) {
  (q * (q - 1) / 4) * log(pi) + sum(lgamma(a + (1 - seq_len(q)) / 2))
}
gwaudit_closed <- function(D, b) {
  q <- nrow(D); a <- (b + q - 1) / 2
  q * (b + q - 1) / 2 * log(2) + gwaudit_log_mgamma_exact(a, q) -
    a * 2 * sum(log(diag(chol(D))))
}

# Independent strict upper-triangular adjacency of an induced subgraph.
gwaudit_upper_adj <- function(g, idx) {
  sub <- igraph::induced_subgraph(g, idx)
  el <- igraph::as_edgelist(sub, names = FALSE); q <- igraph::vcount(sub)
  a <- matrix(0, q, q)
  if (nrow(el)) for (r in seq_len(nrow(el))) {
    a[min(el[r, ]), max(el[r, ])] <- 1
  }
  a
}

# Literal, independent composition over an atom_decomp object.
gwaudit_manual_composition <- function(g, D, b, dec, iter = 100L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  total <- 0
  for (i in seq_along(dec$atoms_index)) {
    idx <- sort(dec$atoms_index[[i]])
    total <- total + if (isTRUE(dec$complete[i])) gwaudit_closed(D[idx, idx, drop = FALSE], b) else
      BDgraph::gnorm(gwaudit_upper_adj(g, idx), b = b, D = D[idx, idx, drop = FALSE],
                     iter = as.integer(iter))
  }
  for (j in seq_along(dec$separators_index)) {
    idx <- sort(dec$separators_index[[j]])
    total <- total - gwaudit_closed(D[idx, idx, drop = FALSE], b)
  }
  total
}

test_that("Path graph: composition equals the literal atom-minus-separator identity", {
  g <- gwaudit_graph("path4"); D <- gwaudit_spd(4, seed = 11)
  fit <- gwaudit_fun("gwishart_log_norm")(g, b = 3, D = D)
  dec <- atom_decompose(g)
  expect_equal(length(dec$atoms_index), 3L)
  expect_equal(length(dec$separators_index), 2L)
  # separators are the single shared vertex of each consecutive edge
  expect_identical(dec$separators_index, list(2L, 3L))
  expect_equal(fit$log_norm, gwaudit_manual_composition(g, D, 3, dec), tolerance = 1e-12)
  expect_false(fit$approximate)
})

test_that("Two overlapping complete cliques: separator counted exactly once with a minus sign", {
  g <- gwaudit_graph("cliques3"); D <- gwaudit_spd(5, seed = 13)
  dec <- atom_decompose(g)
  fit <- gwaudit_fun("gwishart_log_norm")(g, b = 4, D = D)
  expect_true(all(as.logical(dec$complete)))
  expect_identical(dec$atoms_index, list(1:3, 3:5))
  expect_identical(dec$separators_index, list(3L))  # shared vertex only
  manual <- gwaudit_closed(D[1:3, 1:3, drop = FALSE], 4) +
    gwaudit_closed(D[3:5, 3:5, drop = FALSE], 4) -
    gwaudit_closed(D[3, 3, drop = FALSE], 4)
  expect_equal(fit$log_norm, manual, tolerance = 1e-12)
  # independent cross-check through BDgraph's complete branch
  skip_if_not_installed("BDgraph")
  manual_bd <- BDgraph::gnorm(gwaudit_upper_adj(g, 1:3), b = 4, D = D[1:3, 1:3, drop = FALSE], iter = 100L) +
    BDgraph::gnorm(gwaudit_upper_adj(g, 3:5), b = 4, D = D[3:5, 3:5, drop = FALSE], iter = 100L) -
    BDgraph::gnorm(gwaudit_upper_adj(g, 3), b = 4, D = D[3, 3, drop = FALSE], iter = 100L)
  expect_equal(fit$log_norm, manual_bd, tolerance = 1e-10)
})

test_that("A longer chordal chain: every separator is subtracted once, multiplicity respected", {
  g <- gwaudit_graph("path6"); D <- gwaudit_spd(6, seed = 17)
  dec <- atom_decompose(g)
  fit <- gwaudit_fun("gwishart_log_norm")(g, b = 3, D = D)
  expect_identical(length(dec$atoms_index), 5L)
  expect_identical(length(dec$separators_index), 4L)
  expect_equal(fit$log_norm, gwaudit_manual_composition(g, D, 3, dec), tolerance = 1e-12)
  # no separator may appear twice and none may be missing
  keys <- vapply(dec$separators_index, function(s) paste(sort(s), collapse = ","), "")
  expect_false(anyDuplicated(keys) > 0L)
  expect_identical(sort(keys), sort(as.character(2:5)))
  expect_identical(fit$diagnostics$n_separators, length(dec$separators_index))
  expect_equal(sum(fit$separator_log_norm), sum(vapply(fit$separator_fits, function(z) z$log_norm, 0.0)),
               tolerance = 1e-15)
})

test_that("Disconnected complete components: product over components, no artificial separator", {
  g <- gwaudit_graph("disc_complete"); D <- gwaudit_spd(5, seed = 19)
  fit <- gwaudit_fun("gwishart_log_norm")(g, b = 3, D = D)
  expect_identical(fit$diagnostics$n_components, 2L)
  expect_identical(fit$diagnostics$n_separators, 0L)
  expect_equal(length(fit$separator_fits), 0L)
  manual <- gwaudit_closed(D[1:3, 1:3, drop = FALSE], 3) +
    gwaudit_closed(D[4:5, 4:5, drop = FALSE], 3)
  expect_equal(fit$log_norm, manual, tolerance = 1e-12)
  expect_equal(sum(fit$atom_log_norm), manual, tolerance = 1e-12)
})

test_that("Multiple singleton components: each contributes its own closed form", {
  g <- gwaudit_graph("disc_singletons"); D <- gwaudit_spd(3, seed = 23)
  fit <- gwaudit_fun("gwishart_log_norm")(g, b = 5, D = D)
  expect_identical(fit$diagnostics$n_components, 3L)
  expect_identical(fit$diagnostics$n_separators, 0L)
  manual <- sum(vapply(1:3, function(i) gwaudit_closed(D[i, i, drop = FALSE], 5), 0.0))
  expect_equal(fit$log_norm, manual, tolerance = 1e-12)
})

test_that("Disconnected mixed components with aligned RNG", {
  skip_if_not_installed("BDgraph")
  g <- gwaudit_graph("disc_mixed"); D <- gwaudit_spd(7, seed = 29)
  dec <- atom_decompose(g)
  fit <- { set.seed(313); gwaudit_fun("gwishart_log_norm")(g, b = 4, D = D) }
  manual <- gwaudit_manual_composition(g, D, 4, dec, seed = 313)
  expect_identical(fit$diagnostics$n_components, 2L)
  expect_identical(fit$diagnostics$n_separators, 0L)
  expect_true(fit$approximate)
  expect_equal(fit$log_norm, manual, tolerance = 1e-12)
})

test_that("Off-component entries of D do not enter any local computation", {
  g <- gwaudit_graph("disc_complete"); D <- gwaudit_spd(5, seed = 31)
  base <- gwishart_log_norm(g, b = 3, D = D)$log_norm
  # small perturbations of the entries that link the two components; they must
  # leave the value unchanged because only principal blocks are ever used
  Dalt <- D
  eps <- 0.3
  Dalt[1, 4] <- Dalt[4, 1] <- D[1, 4] + eps
  Dalt[1, 5] <- Dalt[5, 1] <- D[1, 5] - eps
  Dalt[3, 5] <- Dalt[5, 3] <- D[3, 5] + eps
  expect_true(all(eigen(Dalt, symmetric = TRUE, only.values = TRUE)$values > 0))
  expect_equal(gwishart_log_norm(g, b = 3, D = Dalt)$log_norm, base, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(Dalt, D, tolerance = 1e-3)))
})

test_that("Component permutation leaves the value invariant", {
  g <- gwaudit_graph("disc_complete"); D <- gwaudit_spd(5, seed = 37)
  base <- gwishart_log_norm(g, b = 3, D = D)$log_norm
  # rebuild the same graph with the components swapped and D permuted accordingly
  g2 <- igraph::set_vertex_attr(
    igraph::disjoint_union(igraph::make_full_graph(2), igraph::make_full_graph(3)),
    "name", value = LETTERS[1:5])
  perm <- c(4, 5, 1, 2, 3)
  D2 <- D[perm, perm]
  expect_equal(gwishart_log_norm(g2, b = 3, D = D2)$log_norm, base, tolerance = 1e-12)
})

test_that("Supplied decomposition is reused and changes nothing numerically", {
  g <- gwaudit_graph("cliques3"); D <- gwaudit_spd(5, seed = 41)
  dec <- atom_decompose(g)
  a <- gwishart_log_norm(g, b = 3, D = D)
  b <- gwishart_log_norm(g, b = 3, D = D, decomp = dec)
  expect_identical(a$log_norm, b$log_norm)
  expect_true(b$diagnostics$decomp_reused)
  expect_false(a$diagnostics$decomp_reused)
  expect_identical(b$timing$decomposition, 0)
  expect_identical(b$decomp$atoms_index, dec$atoms_index)
  expect_identical(b$decomp$separators_index, dec$separators_index)
  # approximate case with aligned RNG
  skip_if_not_installed("BDgraph")
  gc4 <- gwaudit_graph("cycle4"); D4 <- gwaudit_spd(4, seed = 43); dec4 <- atom_decompose(gc4)
  x <- { set.seed(99); gwishart_log_norm(gc4, b = 3, D = D4) }
  y <- { set.seed(99); gwishart_log_norm(gc4, b = 3, D = D4, decomp = dec4) }
  expect_identical(x$log_norm, y$log_norm)
  expect_identical(y$timing$decomposition, 0)
})

test_that("timing$decomposition isolates atom_decompose and is 0 when reused", {
  g <- gwaudit_graph("cliques3"); D <- gwaudit_spd(5, seed = 59)
  dec <- atom_decompose(g)

  # decomp = NULL: only atom_decompose(graph) is timed
  t0 <- proc.time()[[3L]]
  a <- gwishart_log_norm(g, b = 3, D = D)
  expected <- unname(proc.time()[[3L]] - t0)
  expect_false(a$diagnostics$decomp_reused)
  expect_true(is.numeric(a$timing$decomposition))
  expect_length(a$timing$decomposition, 1L)
  expect_true(is.finite(a$timing$decomposition))
  expect_true(a$timing$decomposition >= 0)
  # a single decomposition cannot cost more than the whole call. proc.time()
  # advances in 10 ms steps on Windows, so a bound derived from a separate
  # measurement of the same work may differ by up to one clock tick; the
  # tolerance covers that granularity without weakening the containment claim.
  tick <- 0.02
  expect_lte(a$timing$decomposition, a$timing$total + tick)
  expect_lte(a$timing$decomposition, expected + tick)

  # supplied decomposition: exactly zero, and the strict validation still runs
  b <- gwishart_log_norm(g, b = 3, D = D, decomp = dec)
  expect_identical(b$timing$decomposition, 0)
  expect_true(b$diagnostics$decomp_reused)

  # the timing fix is not allowed to move any statistical result
  expect_identical(a$log_norm, b$log_norm)
  expect_identical(a$atom_log_norm, b$atom_log_norm)
  expect_identical(a$separator_log_norm, b$separator_log_norm)
  expect_identical(a$approximate, b$approximate)
  expect_identical(a$engine, b$engine)
  expect_identical(a$diagnostics$n_atoms, b$diagnostics$n_atoms)
  expect_identical(a$diagnostics$n_separators, b$diagnostics$n_separators)
  expect_identical(a$diagnostics$n_complete_atoms, b$diagnostics$n_complete_atoms)
  expect_identical(a$diagnostics$n_incomplete_atoms, b$diagnostics$n_incomplete_atoms)
  expect_identical(a$diagnostics$n_bdgraph_calls, b$diagnostics$n_bdgraph_calls)
  # local_summary carries the per-locality results plus a wall-clock `elapsed`
  # column. Only the results are comparable run to run; elapsed is a separate
  # measurement of the same work and may differ by one proc.time() clock tick,
  # so it is validated for sanity rather than for equality.
  expect_identical(a$local_summary[setdiff(names(a$local_summary), "elapsed")],
                   b$local_summary[setdiff(names(b$local_summary), "elapsed")])
  expect_identical(names(a$local_summary), names(b$local_summary))
  expect_true(all(is.finite(a$local_summary$elapsed)))
  expect_true(all(a$local_summary$elapsed >= 0))

  # supplied-decomposition validation is still strict: class, structure,
  # vertex order and edge compatibility are all enforced
  expect_error(gwishart_log_norm(g, b = 3, D = D, decomp = list()), "atom_decomp")
  expect_error(gwishart_log_norm(g, b = 3, D = D,
                                 decomp = atom_decompose(gwaudit_graph("path6"))),
               "vertex names|graph edges|cover")
  dec4 <- atom_decompose(gwaudit_graph("cycle4"))
  expect_error(gwishart_log_norm(g, b = 3, D = D, decomp = dec4),
               "vertex names|graph edges|cover")
})

test_that("timing$decomposition stays 0 for an approximate, reused decomposition", {
  skip_if_not_installed("BDgraph")
  gc4 <- gwaudit_graph("cycle4"); D4 <- gwaudit_spd(4, seed = 61)
  dec4 <- atom_decompose(gc4)
  fit <- gwishart_log_norm(gc4, b = 3, D = D4, decomp = dec4)
  expect_identical(fit$timing$decomposition, 0)
  expect_true(fit$diagnostics$decomp_reused)
  expect_true(fit$approximate)
  expect_identical(fit$diagnostics$n_bdgraph_calls, 1L)
  expect_true(is.finite(fit$log_norm))
  expect_true(all(c("decomposition", "local", "composition", "total") %in%
                    names(fit$timing)))
  expect_length(fit$timing, 4L)
})

test_that("A decomposition from another graph is rejected", {
  g <- gwaudit_graph("mixed")
  D <- gwaudit_spd(6, seed = 47)
  wrong <- atom_decompose(gwaudit_graph("path6"))
  expect_error(gwishart_log_norm(g, b = 3, D = D, decomp = wrong))
  wrong2 <- atom_decompose(gwaudit_graph("cycle4"))
  expect_error(gwishart_log_norm(g, b = 3, D = D, decomp = wrong2))
  expect_error(gwishart_log_norm(g, b = 3, D = D, decomp = list()), "atom_decomp")
})
