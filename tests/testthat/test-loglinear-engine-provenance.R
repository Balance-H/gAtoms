# Focused tests for explicit log-linear engine provenance.
#
# The former "auto" engine (prefer gRim, fall back to stats::loglin) has been
# removed: the requested engine is now used as given and a failure is reported
# as an error. These tests pin that provenance contract and the numerical
# regression values. No production logic is mocked except the external engines
# themselves.

llprov_graph <- function(kind) {
  g <- switch(kind,
    path4 = igraph::make_ring(4, circular = FALSE),
    cycle4 = igraph::make_ring(4),
    cycle5 = igraph::make_ring(5))
  igraph::set_vertex_attr(g, "name", value = LETTERS[seq_len(igraph::vcount(g))])
}

llprov_data <- function(g) {
  grid <- do.call(expand.grid,
    stats::setNames(rep(list(c("lo", "hi")), igraph::vcount(g)), igraph::V(g)$name))
  reps <- rep(3L, nrow(grid))
  reps[seq.int(4L, 7L)] <- 4L
  grid[rep(seq_len(nrow(grid)), reps), , drop = FALSE]
}

test_that("The engine argument exposes only the explicit engines with gRim as default", {
  eng <- formals(fit_loglinear_graph)$engine
  expect_identical(eval(eng), c("gRim", "loglin"))
  expect_identical(eval(eng)[1L], "gRim")
  g <- llprov_graph("cycle4"); dat <- llprov_data(g)
  expect_error(fit_loglinear_graph(dat, g, engine = "auto"), "arg")
  expect_error(fit_loglinear_graph(dat, g, engine = "nonsense"), "arg")
})

test_that("An explicit request never falls back to the other engine", {
  g <- llprov_graph("cycle4")
  dat <- llprov_data(g)

  loglin_calls <- 0L
  local_mocked_bindings(capture.output = function(...) stop("simulated gRim failure"),
                        .package = "utils")
  local_mocked_bindings(loglin = function(...) {
    loglin_calls <<- loglin_calls + 1L
    stop("loglin must not be attempted")
  }, .package = "stats")
  e <- tryCatch(fit_loglinear_graph(dat, g, engine = "gRim"), error = identity)
  expect_s3_class(e, "atom_loglinear_error")
  expect_identical(e$engine, "gRim")
  expect_identical(e$requested_engine, "gRim")
  expect_identical(e$attempted_engines, "gRim")
  expect_identical(names(e$engine_errors), "gRim")
  expect_match(conditionMessage(e), "simulated gRim failure", fixed = TRUE)
  expect_identical(loglin_calls, 0L)

  calls <- 0L
  local_mocked_bindings(loglin = function(...) {
    calls <<- calls + 1L
    stop("simulated loglin failure")
  }, .package = "stats")
  e2 <- tryCatch(fit_loglinear_graph(dat, g, engine = "loglin"), error = identity)
  expect_s3_class(e2, "atom_loglinear_error")
  expect_identical(e2$engine, "loglin")
  expect_identical(e2$requested_engine, "loglin")
  expect_identical(e2$attempted_engines, "loglin")
  expect_identical(names(e2$engine_errors), "loglin")
  expect_match(conditionMessage(e2), "simulated loglin failure", fixed = TRUE)
  expect_identical(calls, 1L)
})

test_that("A successful explicit fit records requested, attempted and used engines", {
  g <- llprov_graph("cycle4")
  dat <- llprov_data(g)
  f <- fit_loglinear_graph(dat, g, engine = "loglin")
  z <- f$local_fits[[1]]
  expect_identical(z$requested_engine, "loglin")
  expect_identical(z$attempted_engines, "loglin")
  expect_identical(z$engine, "loglin")
  expect_false(z$fallback)
  expect_true(is.na(z$fallback_reason))
  expect_identical(z$fallback_message, "")
  expect_identical(f$engine, "loglin")
  expect_identical(summary(f)$requested_engine, "loglin")
  expect_identical(summary(f)$fallback_count, 0L)
  expect_length(summary(f)$fallback_reasons, 0L)
  expect_identical(as.character(summary(f)$engine_provenance$attempted_engines), "loglin")
})

test_that("The default engine is gRim", {
  skip_if_not_installed("gRim")
  g <- llprov_graph("cycle4")
  dat <- llprov_data(g)
  f <- fit_loglinear_graph(dat, g)
  z <- f$local_fits[[1]]
  expect_identical(f$engine, "gRim")
  expect_identical(z$requested_engine, "gRim")
  expect_identical(z$attempted_engines, "gRim")
  expect_identical(z$engine, "gRim")
  expect_identical(z$method, "effloglin")
  expect_false(z$fallback)
  # explicit gRim is identical in provenance
  fg <- fit_loglinear_graph(dat, g, engine = "gRim")
  expect_identical(fg$local_fits[[1]]$requested_engine, "gRim")
  expect_identical(fg$local_fits[[1]]$engine, "gRim")
})

test_that("Complete atoms use closed form with no iterative engine requested", {
  g <- llprov_graph("path4")
  dat <- llprov_data(g)
  local_mocked_bindings(loglin = function(...) stop("engine called"), .package = "stats")
  f <- fit_loglinear_graph(dat, g)
  for (z in f$local_fits) {
    expect_identical(z$engine, "closed_form")
    expect_identical(z$requested_engine, "gRim")
    expect_length(z$attempted_engines, 0L)
    expect_false(z$fallback)
    expect_true(is.na(z$fallback_reason))
  }
})

test_that("Likelihood, deviance and fitted values are numerically unchanged", {
  g <- llprov_graph("path4")
  dat <- llprov_data(g)
  f <- fit_loglinear_graph(dat, g, engine = "loglin")
  # frozen reference values from the audited implementation
  expect_equal(as.numeric(logLik(f)), -143.910597794659, tolerance = 1e-9)
  expect_equal(deviance(f), 0.353353196789, tolerance = 1e-9)
  expect_equal(f$likelihood$saturated, -143.733921196264, tolerance = 1e-9)
  expect_equal(sum(fitted(f, type = "full")), 1, tolerance = 1e-12)
  expect_equal(f$likelihood$saturated, sum(f$suffstats$observed_counts *
      log(f$suffstats$observed_counts / f$suffstats$n)), tolerance = 1e-12)
  # no AIC/BIC/df functionality was added
  expect_null(attr(logLik(f), "df"))
  expect_false(any(c("AIC", "BIC") %in% names(summary(f))))
})
