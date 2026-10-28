# Audit-only tests: the precomputed sufficient-margin input pathway of
# fit_loglinear_graph().

llaudit_graph <- function(kind) {
  g <- switch(kind,
    path4 = igraph::make_ring(4, circular = FALSE),
    cycle4 = igraph::make_ring(4),
    mixed = igraph::add_edges(igraph::add_vertices(igraph::make_ring(4), 2),
                              c(3, 5, 4, 5, 5, 6)),
    disco = igraph::disjoint_union(igraph::make_ring(3, circular = FALSE),
                                   igraph::make_empty_graph(1, directed = FALSE)))
  igraph::set_vertex_attr(g, "name", value = LETTERS[seq_len(igraph::vcount(g))])
}

llaudit_data <- function(g) {
  grid <- do.call(expand.grid,
    stats::setNames(rep(list(c("lo", "hi")), igraph::vcount(g)), igraph::V(g)$name))
  reps <- rep(3L, nrow(grid))
  reps[seq.int(4L, 7L)] <- 4L
  grid[rep(seq_len(nrow(grid)), reps), , drop = FALSE]
}

test_that("margin_counts produces the documented loglinear_suffstats object", {
  g <- llaudit_graph("mixed")
  dat <- llaudit_data(g)
  dec <- atom_decompose(g)
  ss <- margin_counts(dat, c(atoms(dec), separators(dec)))
  expect_s3_class(ss, "loglinear_suffstats")
  expect_true(all(c("n", "variables", "levels", "level_counts", "sets", "margins",
                    "observed_counts", "cells") %in% names(ss)))
  expect_identical(ss$variables, igraph::V(g)$name)
  expect_identical(ss$n, nrow(dat))
  expect_identical(ss$level_counts, lengths(ss$levels))
  expect_identical(unname(ss$sets),
                   unname(c(lapply(dec$atoms_index, function(s) igraph::V(g)$name[s]),
                            lapply(dec$separators_index, function(s) igraph::V(g)$name[s]))))
  expect_length(ss$margins, length(ss$sets))
  for (i in seq_along(ss$sets)) {
    expect_identical(dimnames(ss$margins[[i]]), ss$levels[ss$sets[[i]]])
    expect_equal(sum(ss$margins[[i]]), ss$n)
  }
  expect_equal(sum(ss$observed_counts), ss$n)
})

test_that("The sufficient-margin pathway avoids recomputing margins", {
  g <- llaudit_graph("path4")
  dat <- llaudit_data(g)
  dec <- atom_decompose(g)
  ss <- margin_counts(dat, c(atoms(dec), separators(dec)))
  calls <- 0L
  suppressMessages(trace(".loglinear_margins", tracer = function() calls <<- calls + 1L,
                         where = asNamespace("gAtoms"), print = FALSE))
  on.exit(suppressMessages(untrace(".loglinear_margins", where = asNamespace("gAtoms"))), add = TRUE)
  invisible(fit_loglinear_graph(dat, g, decomp = dec, engine = "loglin"))
  raw_calls <- calls
  invisible(fit_loglinear_graph(ss, g, decomp = dec, engine = "loglin"))
  ss_calls <- calls - raw_calls
  expect_equal(raw_calls, 1L)
  expect_equal(ss_calls, 0L)
})

test_that("Raw and sufficient-margin pathways agree on every identifiable summary", {
  for (kind in c("path4", "cycle4", "mixed", "disco")) {
    g <- llaudit_graph(kind)
    dat <- llaudit_data(g)
    dec <- atom_decompose(g)
    ss <- margin_counts(dat, c(atoms(dec), separators(dec)))
    a <- fit_loglinear_graph(dat, g, decomp = dec, engine = "loglin")
    b <- fit_loglinear_graph(ss, g, decomp = dec, engine = "loglin")
    expect_identical(b$suffstats, ss)
    expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)), tolerance = 1e-12)
    expect_equal(deviance(a), deviance(b), tolerance = 1e-12)
    expect_equal(a$likelihood$saturated, b$likelihood$saturated, tolerance = 1e-12)
    expect_equal(fitted(a), fitted(b), tolerance = 1e-12)
    expect_equal(fitted(a, type = "full"), fitted(b, type = "full"), tolerance = 1e-12)
    expect_identical(a$levels, b$levels)
    expect_identical(a$suffstats$n, b$suffstats$n)
    expect_identical(vapply(a$local_fits, function(z) z$engine, ""),
                     vapply(b$local_fits, function(z) z$engine, ""))
    for (i in seq_along(a$local_fits)) {
      expect_equal(a$local_fits[[i]]$fitted_counts, b$local_fits[[i]]$fitted_counts, tolerance = 1e-12)
      expect_identical(a$local_fits[[i]]$generators, b$local_fits[[i]]$generators)
    }
    for (i in seq_along(a$separator_fits)) {
      expect_equal(a$separator_fits[[i]]$fitted_counts, b$separator_fits[[i]]$fitted_counts, tolerance = 1e-12)
    }
  }
})

test_that("Without observed configuration counts deviance and saturated likelihood are NA, not invented", {
  g <- llaudit_graph("mixed")
  dat <- llaudit_data(g)
  dec <- atom_decompose(g)
  ss <- margin_counts(dat, c(atoms(dec), separators(dec)))
  ss$observed_counts <- NULL
  f <- fit_loglinear_graph(ss, g, decomp = dec, engine = "loglin")
  expect_true(is.na(deviance(f)))
  expect_true(is.na(f$likelihood$saturated))
  expect_true(is.finite(as.numeric(logLik(f))))
  expect_true(grepl("unavailable", f$diagnostics$deviance))
  # random observed-configuration counts cannot be recovered from margins alone
  a <- fit_loglinear_graph(dat, g, decomp = dec, engine = "loglin")
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(a)), tolerance = 1e-12)
  expect_equal(fitted(f), fitted(a), tolerance = 1e-12)
})

test_that("Malformed sufficient statistics are rejected with informative errors", {
  g <- llaudit_graph("path4")
  dat <- llaudit_data(g)
  dec <- atom_decompose(g)
  ss <- margin_counts(dat, c(atoms(dec), separators(dec)))
  expect_error(fit_loglinear_graph(margin_counts(dat, list("A")), g), "lack")
  bad <- ss; bad$margins[[1]][1] <- -1
  expect_error(fit_loglinear_graph(bad, g), class = "atom_loglinear_error")
  bad <- ss; bad$n <- bad$n + 1
  expect_error(fit_loglinear_graph(bad, g, decomp = dec), "Invalid sufficient")
  bad <- ss; bad$levels$A <- rev(bad$levels$A)
  expect_error(fit_loglinear_graph(bad, g, decomp = dec), "Invalid sufficient")
  bad <- ss; bad$variables <- rev(bad$variables)
  expect_error(fit_loglinear_graph(bad, g, decomp = dec), "graph variable order")
  # a margin that violates the frozen level metadata
  bad <- ss; dimnames(bad$margins[[1]]) <- NULL
  expect_error(fit_loglinear_graph(bad, g, decomp = dec), class = "atom_loglinear_error")
  # graph incompatibility: margins for a different graph variable set
  other <- margin_counts(dat, list(c("A", "B")))
  expect_error(fit_loglinear_graph(other, g), "lack|Invalid")
})

test_that("The sufficient-margin object is validated even when decomp is recomputed", {
  g <- llaudit_graph("path4")
  dat <- llaudit_data(g)
  dec <- atom_decompose(g)
  ss <- margin_counts(dat, c(atoms(dec), separators(dec)))
  f <- fit_loglinear_graph(ss, g, engine = "loglin")
  expect_s3_class(f, "atom_loglinear_fit")
  expect_true(is.finite(as.numeric(logLik(f))))
  expect_true(is.finite(deviance(f)))
  expect_equal(f$decomp$atoms_index, dec$atoms_index)
  expect_equal(f$decomp$separators_index, dec$separators_index)
  # the suffstats object is reused verbatim, including n and discovered margins
  expect_identical(f$suffstats, ss)
})

test_that("Analysis-of-deviance quantities are available only with raw full configuration counts", {
  g <- llaudit_graph("cycle4")
  dat <- llaudit_data(g)
  dec <- atom_decompose(g)
  ss <- margin_counts(dat, c(atoms(dec), separators(dec)))
  full <- fit_loglinear_graph(dat, g, decomp = dec, engine = "loglin")
  # the saturated reference needs the sparse observed-configuration counts
  expect_true(is.finite(full$likelihood$saturated))
  expect_equal(full$likelihood$deviance, 2 * (full$likelihood$saturated - as.numeric(logLik(full))),
               tolerance = 1e-12)
  expect_null(full$diagnostics$deviance)
  ss$observed_counts <- NULL
  reduced <- fit_loglinear_graph(ss, g, decomp = dec, engine = "loglin")
  expect_true(is.na(reduced$likelihood$deviance))
  expect_true(is.na(reduced$likelihood$saturated))
  expect_false(is.null(reduced$diagnostics$deviance))
})
