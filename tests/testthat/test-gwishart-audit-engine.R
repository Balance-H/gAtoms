# Audit-only tests for the G-Wishart module: BDgraph adapter, integer-b
# contract, input validation and the result-object contract.

gwaudit_ns <- function() asNamespace("gAtoms")

gwaudit_fun <- function(name) {
  if (exists(name, envir = globalenv(), inherits = FALSE)) return(get(name, envir = globalenv()))
  get(name, envir = gwaudit_ns())
}

gwaudit_graph <- function(kind) {
  g <- switch(kind,
    complete3 = igraph::make_full_graph(3),
    path4 = igraph::make_ring(4, circular = FALSE),
    cycle4 = igraph::make_ring(4),
    cycle5 = igraph::make_ring(5),
    mixed = igraph::graph_from_edgelist(rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1),
                                              c(3, 5), c(5, 6)), directed = FALSE),
    disc = igraph::disjoint_union(igraph::make_full_graph(3),
                                  igraph::make_empty_graph(1, directed = FALSE)))
  igraph::set_vertex_attr(g, "name", value = LETTERS[seq_len(igraph::vcount(g))])
}

gwaudit_spd <- function(q, seed = 1L, off = 0.25) {
  set.seed(seed); A <- matrix(off, q, q); diag(A) <- 1; A + diag(q) * 0.7
}

gwaudit_upper_adj <- function(g, idx) {
  sub <- igraph::induced_subgraph(g, idx)
  el <- igraph::as_edgelist(sub, names = FALSE); q <- igraph::vcount(sub)
  a <- matrix(0, q, q)
  if (nrow(el)) for (r in seq_len(nrow(el))) a[min(el[r, ]), max(el[r, ])] <- 1
  a
}

test_that("Single incomplete prime matches a direct BDgraph::gnorm call (cycle4, cycle5)", {
  skip_if_not_installed("BDgraph")
  for (kind in c("cycle4", "cycle5")) for (b in c(3, 4)) {
    g <- gwaudit_graph(kind); q <- igraph::vcount(g)
    D <- gwaudit_spd(q, seed = q + b)
    dec <- atom_decompose(g)
    expect_identical(length(dec$atoms_index), 1L)     # one prime, no separators
    expect_false(isTRUE(dec$complete[1]))
    expect_identical(length(dec$separators_index), 0L)
    set.seed(2024); fit <- gwishart_log_norm(g, b = b, D = D, control = list(iter = 100L))
    set.seed(2024); ref <- BDgraph::gnorm(gwaudit_upper_adj(g, seq_len(q)), b = b, D = D, iter = 100L)
    expect_identical(fit$log_norm, ref)
    expect_identical(fit$diagnostics$n_bdgraph_calls, 1L)
    expect_true(fit$approximate)
    expect_identical(fit$engine, "decomposition+BDgraph")
  }
})

test_that("The adapter passes gnorm exactly adj, b, D_P, iter with no hidden transform", {
  skip_if_not_installed("BDgraph")
  g <- gwaudit_graph("cycle4"); D <- gwaudit_spd(4, seed = 53)
  ns <- gwaudit_ns()
  seen <- list()
  original <- get("gnorm", envir = asNamespace("BDgraph"))
  calls <- 0L
  wrapped <- function(adj, b = 3, D = diag(ncol(adj)), iter = 100) {
    calls <<- calls + 1L
    seen[[calls]] <<- list(adj = adj, b = b, D = D, iter = iter)
    original(adj = adj, b = b, D = D, iter = iter)
  }
  assignInNamespace("gnorm", wrapped, ns = "BDgraph")
  on.exit(assignInNamespace("gnorm", original, ns = "BDgraph"), add = TRUE)
  set.seed(5); fit <- gwishart_log_norm(g, b = 4, D = D, control = list(iter = 37L))
  expect_identical(calls, 1L)                       # no repeated Monte Carlo calls
  a <- seen[[1]]
  expect_identical(a$b, 4)                          # b unchanged after validation
  expect_identical(a$iter, 37L)                     # iter passed unchanged
  expect_identical(a$D, D)                          # same local block, same order
  # adjacency contract
  expect_true(is.matrix(a$adj) && all(dim(a$adj) == c(4, 4)))
  expect_true(all(a$adj %in% c(0, 1)))
  expect_true(all(diag(a$adj) == 0))
  expect_true(all(a$adj[lower.tri(a$adj)] == 0))
  expect_identical(unname(a$adj), gwaudit_upper_adj(g, 1:4))
})

test_that("The stored adjacency for every local fit is the strict upper-triangular form", {
  skip_if_not_installed("BDgraph")
  g <- gwaudit_graph("mixed"); D <- gwaudit_spd(6, seed = 59)
  set.seed(7); fit <- gwishart_log_norm(g, b = 3, D = D)
  expect_gt(sum(vapply(fit$local_fits, function(z) identical(z$engine, "BDgraph"), logical(1))), 0L)
  for (z in fit$local_fits) {
    a <- z$adjacency
    idx <- match(z$vertices, igraph::V(fit$graph)$name)
    expect_true(is.matrix(a) && all(dim(a) == z$size))
    expect_true(all(a %in% c(0, 1)) && all(diag(a) == 0) && all(a[lower.tri(a)] == 0))
    expect_identical(unname(a), gwaudit_upper_adj(fit$graph, idx))
    expect_identical(rownames(a), z$vertices)
  }
})

test_that("Integer-valued b is accepted and non-integer b is rejected before the engine", {
  g_complete <- gwaudit_graph("complete3"); Dc <- gwaudit_spd(3, seed = 61)
  g_cycle <- gwaudit_graph("cycle4"); D4 <- gwaudit_spd(4, seed = 67)
  for (bv in list(3, 3L, 3.0, 4, 10)) {
    expect_silent(fit <- gwishart_log_norm(g_complete, b = bv, D = Dc))
    expect_identical(fit$b, as.numeric(bv))
    expect_silent(fit2 <- gwishart_log_norm(g_cycle, b = bv, D = D4))
    expect_identical(fit2$b, as.numeric(bv))
  }
  for (bv in list(2, 2.999, 3.1, 3.5, 4.25, NA_real_, NaN, Inf, -Inf,
                  numeric(0), c(3, 4), "3", .Machine$integer.max + 1)) {
    expect_error(gwishart_log_norm(g_complete, b = bv, D = Dc))
    expect_error(gwishart_log_norm(g_cycle, b = bv, D = D4))
  }
  # graph-independent: the complete graph must reject 3.5 exactly as the cycle does
  expect_error(gwishart_log_norm(g_complete, b = 3.5, D = Dc), "integer-valued")
  expect_error(gwishart_log_norm(g_cycle, b = 3.5, D = D4), "integer-valued")
})

test_that("Invalid b is rejected before any local engine work", {
  ns <- gwaudit_ns()
  g <- gwaudit_graph("cycle4"); D <- gwaudit_spd(4, seed = 71)
  original <- get(".gwishart_bdgraph_log_norm", envir = ns)
  assignInNamespace(".gwishart_bdgraph_log_norm",
                    function(...) stop("SENTINEL engine reached"), ns = "gAtoms")
  on.exit(assignInNamespace(".gwishart_bdgraph_log_norm", original, ns = "gAtoms"), add = TRUE)
  expect_error(gwishart_log_norm(g, b = 3.5, D = D), "integer-valued")
  # a valid b does reach the engine, which proves the sentinel is active
  expect_error(gwishart_log_norm(g, b = 3, D = D), "SENTINEL")
})

test_that("The b validator itself enforces the integer contract", {
  vb <- get(".gwishart_validate_b", envir = gwaudit_ns())
  expect_identical(vb(3), 3)
  expect_identical(vb(3L), 3)
  expect_identical(vb(3.0), 3)
  expect_identical(vb(2147483647), 2147483647)
  for (bv in list(2, 2.999, 3.5, NA_real_, NaN, Inf, numeric(0), c(3, 4), "3"))
    expect_error(vb(bv))
  expect_error(vb(.Machine$integer.max + 1), "integer range")
})

test_that("Exact and approximate labelling is consistent", {
  skip_if_not_installed("BDgraph")
  exact <- gwishart_log_norm(gwaudit_graph("path4"), b = 3, D = gwaudit_spd(4, seed = 73))
  expect_false(exact$approximate)
  expect_identical(exact$engine, "closed_form")
  expect_identical(exact$diagnostics$n_bdgraph_calls, 0L)
  expect_true(grepl("all prime components and separators evaluated in closed form",
                    exact$diagnostics$note, fixed = TRUE))
  set.seed(3)
  approx <- gwishart_log_norm(gwaudit_graph("cycle4"), b = 3, D = gwaudit_spd(4, seed = 79))
  expect_true(approx$approximate)
  expect_identical(approx$engine, "decomposition+BDgraph")
  expect_gt(approx$diagnostics$n_bdgraph_calls, 0L)
  expect_true(grepl("approximated by BDgraph Monte Carlo", approx$diagnostics$note, fixed = TRUE))
  # only incomplete primes are flagged approximate
  expect_true(all(vapply(approx$separator_fits, function(z) !isTRUE(z$approximate), logical(1))))
  expect_true(all(vapply(approx$separator_fits, function(z) z$engine == "closed_form", logical(1))))
  # no invented uncertainty fields
  expect_false(any(c("se", "sd", "stderr", "std_error", "variance", "uncertainty",
                     "ci", "confint") %in% names(approx)))
  expect_false(any(grepl("stderr|std_error|variance|uncertainty", names(approx),
                         ignore.case = TRUE)))
})

test_that("Result object fields and counts are internally consistent", {
  skip_if_not_installed("BDgraph")
  g <- gwaudit_graph("mixed"); D <- gwaudit_spd(6, seed = 83)
  set.seed(11); fit <- gwishart_log_norm(g, b = 4, D = D)
  expect_s3_class(fit, "atom_gwishart_norm")
  expect_true(all(c("log_norm", "b", "D", "graph", "decomp", "atom_log_norm",
                    "separator_log_norm", "local_fits", "separator_fits",
                    "local_summary", "approximate", "diagnostics", "timing",
                    "control", "engine", "call") %in% names(fit)))
  expect_equal(fit$log_norm, sum(fit$atom_log_norm) - sum(fit$separator_log_norm),
               tolerance = 1e-15)
  expect_identical(fit$diagnostics$n_atoms, length(fit$local_fits))
  expect_identical(fit$diagnostics$n_separators, length(fit$separator_fits))
  expect_identical(fit$diagnostics$n_complete_atoms,
                   sum(vapply(fit$local_fits, function(z) isTRUE(z$complete), logical(1))))
  expect_identical(fit$diagnostics$n_incomplete_atoms,
                   sum(vapply(fit$local_fits, function(z) !isTRUE(z$complete), logical(1))))
  expect_identical(fit$diagnostics$n_bdgraph_calls,
                   sum(vapply(fit$local_fits, function(z) identical(z$engine, "BDgraph"), logical(1))))
  expect_identical(fit$diagnostics$p, igraph::vcount(fit$graph))
  expect_identical(fit$diagnostics$full_graph_size, igraph::vcount(fit$graph))
  expect_identical(fit$b, 4)
  expect_identical(fit$control$iter, 100L)
  expect_identical(dim(fit$D), c(6L, 6L))
  expect_true(all(vapply(fit$atom_log_norm, is.finite, logical(1))))
  expect_true(all(is.finite(fit$separator_log_norm)))
  expect_identical(nrow(fit$local_summary), length(fit$local_fits) + length(fit$separator_fits))
  expect_true(all(c("decomposition", "local", "composition", "total") %in% names(fit$timing)))
})

test_that("D validation rejects invalid matrices without repair", {
  g <- gwaudit_graph("cycle4"); D <- gwaudit_spd(4, seed = 89)
  expect_error(gwishart_log_norm(g, b = 3, D = as.numeric(D)), "numeric matrix")
  expect_error(gwishart_log_norm(g, b = 3, D = matrix(1, 3, 4)), "square")
  expect_error(gwishart_log_norm(g, b = 3, D = gwaudit_spd(3)), "dimension")
  expect_error(gwishart_log_norm(g, b = 3, D = gwaudit_spd(5)), "dimension")
  Asym <- D; Asym[1, 2] <- Asym[1, 2] + 0.5
  expect_error(gwishart_log_norm(g, b = 3, D = Asym), "symmetric")
  expect_error(gwishart_log_norm(g, b = 3, D = matrix(1, 4, 4)), "positive definite")
  expect_error(gwishart_log_norm(g, b = 3, D = diag(c(1, 1, 1, -2))), "positive definite")
  Dna <- D; Dna[1, 1] <- NA
  expect_error(gwishart_log_norm(g, b = 3, D = Dna), "NA, NaN or infinite")
  Dnan <- D; Dnan[2, 2] <- NaN
  expect_error(gwishart_log_norm(g, b = 3, D = Dnan), "NA, NaN or infinite")
  Dinf <- D; Dinf[3, 3] <- Inf
  expect_error(gwishart_log_norm(g, b = 3, D = Dinf), "NA, NaN or infinite")
  Dninf <- D; Dninf[4, 4] <- -Inf
  expect_error(gwishart_log_norm(g, b = 3, D = Dninf), "NA, NaN or infinite")
})

test_that("The documented symmetry tolerance is 1e-8 relative and does not repair large asymmetry", {
  g <- gwaudit_graph("cycle4"); D <- gwaudit_spd(4, seed = 97)
  tiny <- D; tiny[1, 2] <- tiny[1, 2] + 1e-10 * max(abs(D))
  expect_silent(f0 <- gwishart_log_norm(g, b = 3, D = tiny))
  # the internally symmetrized matrix is used, so the result equals the symmetric part
  Sym <- (tiny + t(tiny)) / 2
  expect_equal(f0$D, Sym, tolerance = 1e-15)
  big <- D; big[1, 2] <- big[1, 2] + 1e-4 * max(abs(D))
  expect_error(gwishart_log_norm(g, b = 3, D = big), "symmetric")
})

test_that("control validation is strict", {
  g <- gwaudit_graph("complete3"); D <- gwaudit_spd(3, seed = 101)
  expect_identical(gwishart_log_norm(g, b = 3, D = D)$control$iter, 100L)
  expect_identical(gwishart_log_norm(g, b = 3, D = D, control = list(iter = 7L))$control$iter, 7L)
  expect_identical(gwishart_log_norm(g, b = 3, D = D, control = list(iter = 7))$control$iter, 7L)
  expect_error(gwishart_log_norm(g, b = 3, D = D, control = list(iter = 0)), "iter")
  expect_error(gwishart_log_norm(g, b = 3, D = D, control = list(iter = -1)), "iter")
  expect_error(gwishart_log_norm(g, b = 3, D = D, control = list(iter = 1.5)), "iter")
  expect_error(gwishart_log_norm(g, b = 3, D = D, control = list(iter = Inf)), "iter")
  expect_error(gwishart_log_norm(g, b = 3, D = D, control = list(unknown = 1)), "recognized options")
  expect_error(gwishart_log_norm(g, b = 3, D = D, control = "x"), "named list")
})

test_that("log_norm increases with b on a fixed graph and D (G-Wishart sanity)", {
  g <- gwaudit_graph("complete3"); D <- gwaudit_spd(3, seed = 103)
  vals <- vapply(c(3, 4, 5, 6), function(bb) gwishart_log_norm(g, b = bb, D = D)$log_norm, 0.0)
  expect_true(all(diff(vals) > 0))
})

test_that("engine argument is validated without altering behaviour", {
  g <- gwaudit_graph("path4"); D <- gwaudit_spd(4, seed = 107)
  a <- gwishart_log_norm(g, b = 3, D = D)
  b <- gwishart_log_norm(g, b = 3, D = D, engine = "BDgraph")
  expect_identical(a$log_norm, b$log_norm)
  expect_error(gwishart_log_norm(g, b = 3, D = D, engine = "closed_form"))
  expect_error(gwishart_log_norm(g, b = 3, D = D, engine = "nonsense"))
})
