# Audit-only tests: likelihood convention of the graphical log-linear module.
# These tests are strictly read-only with respect to the package: they add no
# production code and change no existing file. Every reference quantity is
# computed independently of the gAtoms likelihood helpers (.loglinear_logprob,
# .loglinear_count_loglik) so that validation is not circular.

llaudit_graph <- function(kind) {
  g <- switch(kind,
    path3 = igraph::make_ring(3, circular = FALSE),
    path4 = igraph::make_ring(4, circular = FALSE),
    cycle4 = igraph::make_ring(4),
    cycle5 = igraph::make_ring(5),
    complete3 = igraph::make_full_graph(3),
    mixed = igraph::add_edges(igraph::add_vertices(igraph::make_ring(4), 2),
                              c(3, 5, 4, 5, 5, 6)),
    disco = igraph::disjoint_union(igraph::make_ring(3, circular = FALSE),
                                   igraph::make_empty_graph(1, directed = FALSE)))
  igraph::set_vertex_attr(g, "name", value = LETTERS[seq_len(igraph::vcount(g))])
}

# Deterministically supported data: consecutive blocks of the full Cartesian grid
# are repeated so that every cell of every generator margin is observed.
llaudit_data <- function(g, base = 3L, rep_from = 4L, rep_to = 7L) {
  grid <- do.call(expand.grid,
    stats::setNames(rep(list(c("lo", "hi")), igraph::vcount(g)), igraph::V(g)$name))
  reps <- rep(base, nrow(grid))
  reps[seq.int(rep_from, rep_to)] <- base + 1L
  grid[rep(seq_len(nrow(grid)), reps), , drop = FALSE]
}

# Independent reference: fit the model with generators = maximal cliques of G
# using stats::loglin on a genuinely materialized full table, then evaluate the
# multinomial probability log-likelihood DIRECTLY from the returned fitted
# counts. No gAtoms helper participates in this calculation.
llaudit_reference <- function(dat, g, levels) {
  tab <- do.call(table, Map(function(x, l) factor(as.character(x), levels = l), dat, levels))
  n <- sum(tab)
  counts <- as.vector(tab)
  positive <- counts > 0
  z <- stats::loglin(tab, margin = lapply(igraph::max_cliques(g), as.integer),
                     fit = TRUE, eps = 1e-12, iter = 20000, print = FALSE)
  fitted_counts <- as.vector(z$fit)
  ll <- sum(counts[positive] * log(fitted_counts[positive] / n))
  sat <- sum(counts[positive] * log(counts[positive] / n))
  list(n = n, counts = counts, fitted_counts = fitted_counts, ll = ll, sat = sat,
       deviance = 2 * (sat - ll), prob = fitted_counts / n)
}

test_that("logLik is the multinomial probability log-likelihood without factorial constants", {
  for (kind in c("path3", "path4", "cycle4", "mixed", "disco")) {
    g <- llaudit_graph(kind)
    dat <- llaudit_data(g)
    f <- fit_loglinear_graph(dat, g, engine = "loglin")
    ref <- llaudit_reference(dat, g, f$levels)
    # absolute logLik agrees with the independently computed reference
    expect_equal(as.numeric(logLik(f)), ref$ll, tolerance = 1e-8)
    # It is NOT the raw Poisson log-likelihood sum_x n_x log m_x - m_x - log(n_x!).
    # With m_x = n p_x (so sum_x m_x = n) the exact, model-independent identity is
    #   sum_x n_x log m_x - m_x - log(n_x!)
    #     = sum_x n_x log p_x - n log(n) + n + sum_x log(n_x!)
    # i.e. l_mult = l_pois - n*log(n) + n + sum_x log(n_x!).
    poisson_ll <- sum(ref$counts * log(ref$fitted_counts) - ref$fitted_counts -
                        lgamma(ref$counts + 1))
    expect_equal(as.numeric(logLik(f)),
                 poisson_ll - ref$n * log(ref$n) + ref$n + sum(lgamma(ref$counts + 1)),
                 tolerance = 1e-8)
    # The full multinomial log-likelihood including its combinatorial constant is
    # then l_mult + sum_x log(n_x!) = l_pois - n log(n) + n + 2 sum_x log(n_x!).
    expect_equal(ref$ll + sum(lgamma(ref$counts + 1)),
                 poisson_ll - ref$n * log(ref$n) + ref$n + 2 * sum(lgamma(ref$counts + 1)),
                 tolerance = 1e-8)
  }
})

test_that("Saturated likelihood of a complete atom reproduces the empirical entropy", {
  g <- llaudit_graph("complete3")
  dat <- llaudit_data(g)
  tab <- table(dat)
  counts <- as.vector(tab)
  n <- sum(tab)
  positive <- counts > 0
  empirical <- sum(counts[positive] * log(counts[positive] / n))
  f <- fit_loglinear_graph(dat, g, engine = "loglin")
  expect_true(all(f$decomp$complete))
  expect_equal(as.numeric(logLik(f)), empirical, tolerance = 1e-10)
  expect_equal(f$likelihood$saturated, empirical, tolerance = 1e-10)
  expect_equal(deviance(f), 0)
  expect_equal(as.vector(fitted(f, type = "full")), as.vector(counts / n), tolerance = 1e-12)
})

test_that("Deviance equals the likelihood-ratio statistic of an external implementation", {
  for (kind in c("path4", "cycle4", "cycle5", "mixed", "disco")) {
    g <- llaudit_graph(kind)
    dat <- llaudit_data(g)
    f <- fit_loglinear_graph(dat, g, engine = "loglin")
    ref <- llaudit_reference(dat, g, f$levels)
    expect_equal(deviance(f), ref$deviance, tolerance = 1e-8)
    # internal consistency: deviance = 2*(saturated - logLik) under ONE convention
    expect_equal(deviance(f), 2 * (f$likelihood$saturated - as.numeric(logLik(f))),
                 tolerance = 1e-10)
    # independent direct reconstruction from the materialized fitted table
    tab <- table(dat)
    counts <- as.vector(tab)
    positive <- counts > 0
    fitprob <- as.vector(fitted(f, type = "full"))
    expect_equal(deviance(f), 2 * sum(counts[positive] * log(counts[positive] /
      (ref$n * fitprob[positive]))), tolerance = 1e-8)
  }
})

test_that("MASS::loglm agrees on deviance and the multinomial logLik, not on the Poisson constant", {
  skip_if_not_installed("MASS")
  g <- llaudit_graph("mixed")
  dat <- llaudit_data(g)
  f <- fit_loglinear_graph(dat, g, engine = "loglin")
  tab <- table(dat)
  counts <- as.vector(tab)
  n <- sum(tab)
  positive <- counts > 0
  # generating class of the mixed graph (its maximal cliques)
  oracle <- MASS::loglm(~ A * B + A * D + B * C + C * D * E + E * F,
                        data = tab, fit = TRUE, eps = 1e-12, iter = 20000)
  oracle_counts <- as.vector(oracle$fit)
  oracle_ll <- sum(counts[positive] * log(oracle_counts[positive] / n))
  expect_equal(deviance(f), as.numeric(oracle$lrt), tolerance = 1e-7)
  expect_equal(as.numeric(logLik(f)), oracle_ll, tolerance = 1e-7)
  # loglm reports the Poisson deviance; the two log-likelihoods differ exactly by
  # the model-independent constant
  expect_equal(oracle_ll - as.numeric(logLik(f)), 0, tolerance = 1e-7)
  expect_equal(as.vector(fitted(f, type = "full")), oracle_counts / n, tolerance = 1e-7)
})

test_that("logLik carries class logLik and nobs, and the fitted object exposes no df/AIC", {
  g <- llaudit_graph("path4")
  dat <- llaudit_data(g)
  f <- fit_loglinear_graph(dat, g, engine = "loglin")
  ll <- logLik(f)
  expect_s3_class(ll, "logLik")
  expect_identical(attr(ll, "nobs"), f$suffstats$n)
  expect_identical(as.numeric(ll), f$likelihood$logLik[[1]])
  # documented limitation: model df and information criteria are not provided
  expect_null(attr(ll, "df"))
  s <- summary(f)
  expect_false(any(c("AIC", "BIC", "df") %in% names(s)))
  expect_true("logLik" %in% names(s))
  expect_true("deviance" %in% names(s))
})

test_that("Atom and separator contributions use one convention and reproduce global logLik", {
  g <- llaudit_graph("mixed")
  dat <- llaudit_data(g)
  f <- fit_loglinear_graph(dat, g, engine = "loglin")
  # recompute the factorization by hand from stored probabilities (independent of
  # .loglinear_count_loglik)
  term <- function(z) {
    counts <- as.vector(z$observed)
    prob <- as.vector(z$fitted_prob)
    positive <- counts > 0
    sum(counts[positive] * log(prob[positive]))
  }
  manual <- sum(vapply(f$local_fits, term, 0)) - sum(vapply(f$separator_fits, term, 0))
  expect_equal(manual, as.numeric(logLik(f)), tolerance = 1e-10)
  # the saturated term uses the same probability convention
  expect_equal(f$likelihood$saturated,
               sum(as.vector(f$suffstats$observed_counts) *
                     log(as.vector(f$suffstats$observed_counts) / f$suffstats$n)),
               tolerance = 1e-10)
})

test_that("Local fits reproduce observed generator margins exactly (IPF contract)", {
  for (kind in c("cycle4", "cycle5", "mixed")) {
    g <- llaudit_graph(kind)
    dat <- llaudit_data(g)
    f <- fit_loglinear_graph(dat, g, engine = "loglin")
    for (z in f$local_fits) {
      for (gen in z$generators) {
        s <- match(gen, z$vertices)
        expect_equal(as.vector(margin.table(z$fitted_prob, s)),
                     as.vector(margin.table(z$observed, s)) / sum(z$observed),
                     tolerance = 1e-6)
      }
      expect_equal(sum(z$fitted_prob), 1, tolerance = 1e-8)
    }
    for (z in f$separator_fits) expect_equal(sum(z$fitted_prob), 1, tolerance = 1e-8)
  }
})
