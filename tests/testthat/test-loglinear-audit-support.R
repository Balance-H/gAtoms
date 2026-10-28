# Audit-only tests: zero separator support and zero/zero configurations in the
# factorized graphical log-linear representation.

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

llaudit_data <- function(g) {
  grid <- do.call(expand.grid,
    stats::setNames(rep(list(c("lo", "hi")), igraph::vcount(g)), igraph::V(g)$name))
  reps <- rep(3L, nrow(grid))
  reps[seq.int(4L, 7L)] <- 4L
  grid[rep(seq_len(nrow(grid)), reps), , drop = FALSE]
}

llaudit_logprob <- function(fac, data = NULL, codes = NULL) {
  get(".loglinear_logprob", asNamespace("gAtoms"))(fac, data = data, codes = codes)
}

# Independent re-implementation of the factorized log-probability, used to cross
# check the production evaluator on configurations outside separator support.
llaudit_manual_logprob <- function(fac, grid) {
  m <- nrow(grid)
  out <- numeric(m)
  zero <- rep(FALSE, m)
  push <- function(vertices, prob, sign) {
    s <- match(vertices, fac$variables)
    idx <- rep(1L, m)
    stride <- c(1, cumprod(as.double(dim(prob))))[seq_along(dim(prob))]
    for (j in seq_along(s)) {
      code <- match(as.character(grid[[s[j]]]), dimnames(prob)[[j]])
      idx <- idx + (code - 1L) * stride[j]
    }
    p <- as.vector(prob)[idx]
    zero <<- zero | p == 0
    pos <- p > 0
    out[pos] <<- out[pos] + sign * log(p[pos])
  }
  for (z in fac$atoms) push(z$vertices, z$prob, 1)
  for (z in fac$separators) push(z$vertices, z$prob, -1)
  out[zero] <- -Inf
  out
}

test_that("A zero separator marginal yields probability zero, never NaN or Inf", {
  g <- llaudit_graph("path4")
  dat <- llaudit_data(g)
  dat <- dat[dat$B == "hi", , drop = FALSE]          # B = "lo" never observed
  expect_gt(nrow(dat), 0L)
  f <- fit_loglinear_graph(dat, g, engine = "loglin")
  full <- fitted(f, type = "full")
  expect_false(anyNA(full))
  expect_true(all(is.finite(full)))
  expect_true(all(full >= 0))
  expect_equal(sum(full), 1, tolerance = 1e-10)
  # the unsupported separator state carries no mass
  expect_equal(sum(full[, "lo", , ]), 0, tolerance = 1e-12)
  # the separator factor really has a zero entry
  expect_true(any(vapply(f$separator_fits, function(z) any(z$fitted_prob == 0), logical(1))))
  # factorized log-probability stays valid on the whole Cartesian grid
  grid <- do.call(expand.grid, f$levels)
  lp <- llaudit_logprob(fitted(f), grid)
  expect_false(any(is.nan(lp)))
  expect_false(any(lp == Inf))
  expect_true(all(is.finite(lp) | lp == -Inf))
  expect_equal(sum(exp(lp)), 1, tolerance = 1e-10)
  expect_equal(lp, llaudit_manual_logprob(fitted(f), grid), tolerance = 1e-12)
})

test_that("Zero-over-zero is structurally impossible: separator and atom zero supports agree", {
  g <- llaudit_graph("path4")
  dat <- llaudit_data(g)
  dat <- dat[dat$B == "hi", , drop = FALSE]
  f <- fit_loglinear_graph(dat, g, engine = "loglin")
  checked <- 0L
  for (sep in f$separator_fits) {
    for (z in f$local_fits) {
      s <- match(sep$vertices, z$vertices)
      if (anyNA(s)) next
      atom_margin <- margin.table(z$fitted_prob, s)
      expect_identical(dimnames(atom_margin), dimnames(sep$fitted_prob))
      # whenever the separator has zero mass the neighbouring atom does too, and
      # vice versa; hence -Inf - (-Inf) cannot arise through this pair
      expect_identical(atom_margin == 0, sep$fitted_prob == 0)
      checked <- checked + 1L
    }
  }
  expect_gt(checked, 0L)
})

test_that("Factorized and full output agree on zero-support configurations", {
  g <- llaudit_graph("path4")
  dat <- llaudit_data(g)
  dat <- dat[dat$B == "hi", , drop = FALSE]
  f <- fit_loglinear_graph(dat, g, engine = "loglin")
  full <- fitted(f, type = "full")
  grid <- do.call(expand.grid, f$levels)
  expect_equal(as.vector(full), exp(llaudit_logprob(fitted(f), grid)), tolerance = 1e-12)
  zero_cells <- as.vector(full) == 0
  expect_gt(sum(zero_cells), 0L)
  lp <- llaudit_logprob(fitted(f), grid)
  expect_true(all(lp[zero_cells] == -Inf))
})

test_that("Declared unused levels are preserved and receive zero probability", {
  g <- llaudit_graph("path4")
  dat <- llaudit_data(g)
  dat$A <- factor(dat$A, levels = c("unused", "lo", "hi"))
  dat <- dat[dat$A != "unused", , drop = FALSE]
  f <- fit_loglinear_graph(dat, g, engine = "loglin")
  expect_identical(f$levels$A, c("unused", "lo", "hi"))
  full <- fitted(f, type = "full")
  expect_equal(sum(full["unused", , , ]), 0, tolerance = 1e-12)
  expect_false(anyNA(full))
  expect_equal(sum(full), 1, tolerance = 1e-10)
  # unused levels are still in every margin's metadata, with zero counts
  for (z in c(f$local_fits, f$separator_fits)) {
    if ("A" %in% z$vertices) {
      s <- match("A", z$vertices)
      expect_equal(sum(margin.table(z$observed, s)[dimnames(z$observed)[[s]] == "unused"]), 0)
    }
  }
})

test_that("Sampling zeros are not confused with structural or unused-level zeros", {
  g <- llaudit_graph("path4")
  dat <- llaudit_data(g)
  count_before <- nrow(dat) - sum(dat$A == "lo" & dat$B == "lo")
  dat <- dat[!(dat$A == "lo" & dat$B == "lo"), , drop = FALSE]
  expect_equal(nrow(dat), count_before)
  f <- fit_loglinear_graph(dat, g, engine = "loglin")
  full <- fitted(f, type = "full")
  # a sampling zero inside a fully supported atom stays at probability zero here
  expect_equal(sum(full["lo", "lo", , ]), 0, tolerance = 1e-12)
  expect_false(anyNA(full))
  expect_equal(sum(full), 1, tolerance = 1e-10)
  # no structural-zero model is claimed: the fit itself is a plain saturated
  # closed form for a complete-atom graph and simply reproduces the table
  g2 <- llaudit_graph("complete3")
  dat2 <- llaudit_data(g2)
  dat2 <- dat2[!(dat2$A == "lo" & dat2$B == "lo" & dat2$C == "lo"), , drop = FALSE]
  f2 <- fit_loglinear_graph(dat2, g2, engine = "loglin")
  expect_equal(sum(fitted(f2, type = "full")["lo", "lo", "lo"]), 0, tolerance = 1e-12)
  expect_equal(deviance(f2), 0, tolerance = 1e-10)
})

test_that("Observed positive configurations always keep positive fitted probability", {
  for (kind in c("path4", "cycle4", "mixed")) {
    g <- llaudit_graph(kind)
    dat <- llaudit_data(g)
    dat <- dat[stats::runif(nrow(dat)) < 0.6, , drop = FALSE]
    f <- fit_loglinear_graph(dat, g, engine = "loglin")
    lp <- llaudit_logprob(fitted(f), dat)
    expect_true(all(is.finite(lp)))
    expect_true(all(lp < 0))
    # logLik is finite and equals the observed-configuration sum
    expect_true(is.finite(as.numeric(logLik(f))))
  }
})

test_that("Evaluation of the factorized object is robust to inconsistent factors", {
  g <- llaudit_graph("path4")
  dat <- llaudit_data(g)
  dat <- dat[dat$B == "hi", , drop = FALSE]
  f <- fit_loglinear_graph(dat, g, engine = "loglin")
  fac <- fitted(f)
  grid <- do.call(expand.grid, fac$levels)
  i <- which(vapply(fac$separators, function(z) "B" %in% z$vertices, logical(1)))[1]
  # force an all-zero separator factor: every configuration must become -Inf,
  # and 0/0 must not leak through as NaN
  fac0 <- fac
  fac0$separators[[i]]$prob[] <- 0
  lp0 <- llaudit_logprob(fac0, grid)
  expect_false(any(is.nan(lp0)))
  expect_true(all(lp0 == -Inf))
  # force a partly-zero separator factor while atoms stay positive: the zero
  # states become -Inf and no NaN appears
  fac1 <- fac
  fac1$separators[[i]]$prob[] <- rep(c(0, 1), length.out = length(fac1$separators[[i]]$prob))
  lp1 <- llaudit_logprob(fac1, grid)
  expect_false(any(is.nan(lp1)))
  expect_true(any(lp1 == -Inf))
  expect_true(any(is.finite(lp1)))
})

test_that("Full-table materialization stays normalized under extreme sparsity and disconnection", {
  cases <- list(
    constant_B = local({
      g <- llaudit_graph("path4")
      list(g = g, d = data.frame(A = factor(rep(c("lo", "hi"), 20)), B = factor(rep("hi", 40)),
                                 C = factor(rep(c("u", "v"), 20)), D = factor(rep(c("p", "q"), 20))))
    }),
    single_row = local({
      g <- llaudit_graph("path4")
      list(g = g, d = data.frame(A = factor("lo"), B = factor("hi"), C = factor("u"), D = factor("p")))
    }),
    three_rows = local({
      g <- llaudit_graph("disco")
      list(g = g, d = data.frame(A = factor(c("lo", "hi", "lo")), B = factor(c("hi", "hi", "lo")),
                                 C = factor(c("u", "v", "u")), D = factor(c("p", "q", "p"))))
    })
  )
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    f <- fit_loglinear_graph(cs$d, cs$g, engine = "loglin")
    full <- fitted(f, type = "full")
    expect_false(anyNA(full))
    expect_true(all(is.finite(full)))
    expect_true(all(full >= 0))
    expect_equal(sum(full), 1, tolerance = 1e-8)
    lp <- llaudit_logprob(fitted(f), cs$d)
    expect_true(all(is.finite(lp)))
  }
})
