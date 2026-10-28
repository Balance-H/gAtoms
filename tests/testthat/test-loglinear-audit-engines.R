# Audit-only tests: explicit local-engine selection for the log-linear module.
# The former "auto" engine (prefer gRim, fall back to stats::loglin) has been
# removed, so these tests cover only the explicit engines "gRim" and "loglin".
# They mock only external dependencies (stats::loglin, utils::capture.output);
# no production logic is modified.

llaudit_graph <- function(kind) {
  g <- switch(kind,
    path4 = igraph::make_ring(4, circular = FALSE),
    cycle4 = igraph::make_ring(4),
    cycle5 = igraph::make_ring(5),
    mixed = igraph::add_edges(igraph::add_vertices(igraph::make_ring(4), 2),
                              c(3, 5, 4, 5, 5, 6)))
  igraph::set_vertex_attr(g, "name", value = LETTERS[seq_len(igraph::vcount(g))])
}

llaudit_data <- function(g) {
  grid <- do.call(expand.grid,
    stats::setNames(rep(list(c("lo", "hi")), igraph::vcount(g)), igraph::V(g)$name))
  reps <- rep(3L, nrow(grid))
  reps[seq.int(4L, 7L)] <- 4L
  grid[rep(seq_len(nrow(grid)), reps), , drop = FALSE]
}

test_that("Complete atoms always use the closed form and never call an engine", {
  g <- llaudit_graph("path4")
  dat <- llaudit_data(g)
  expect_true(all(atom_decompose(g)$complete))
  local_mocked_bindings(loglin = function(...) stop("engine called"),
                        .package = "stats")
  f <- fit_loglinear_graph(dat, g)
  expect_true(all(vapply(f$local_fits, function(z) z$engine == "closed_form", logical(1))))
  expect_true(all(vapply(f$local_fits, function(z) z$method == "saturated", logical(1))))
  expect_identical(vapply(f$local_fits, function(z) z$iterations, 0L), rep(0L, length(f$local_fits)))
  # separators are always closed-form saturated fits
  expect_true(all(vapply(f$separator_fits, function(z) z$engine == "closed_form", logical(1))))
})

test_that("engine='gRim' uses effloglin and agrees with engine='loglin'", {
  skip_if_not_installed("gRim")
  g <- llaudit_graph("cycle4")
  dat <- llaudit_data(g)
  f <- fit_loglinear_graph(dat, g, engine = "gRim")
  expect_false(all(atom_decompose(g)$complete))
  expect_true(any(vapply(f$local_fits, function(z) z$engine == "gRim", logical(1))))
  expect_true(any(vapply(f$local_fits, function(z) z$method == "effloglin", logical(1))))
  expect_identical(f$engine, "gRim")
  expect_true(all(vapply(f$local_fits, function(z) z$requested_engine == "gRim", logical(1))))
  expect_true(all(vapply(f$local_fits, function(z) identical(z$attempted_engines, "gRim"), logical(1))))
  expect_true(all(vapply(f$local_fits, function(z) !isTRUE(z$fallback), logical(1))))
  # gRim and loglin must agree on the fitted probabilities
  fl <- fit_loglinear_graph(dat, g, engine = "loglin")
  for (i in seq_along(f$local_fits)) {
    expect_equal(f$local_fits[[i]]$fitted_prob, fl$local_fits[[i]]$fitted_prob, tolerance = 1e-6)
  }
})

test_that("engine='loglin' only uses stats::loglin, with no fallback to gRim", {
  g <- llaudit_graph("cycle4")
  dat <- llaudit_data(g)
  grim_calls <- 0L
  local_mocked_bindings(effloglin = function(...) {
    grim_calls <<- grim_calls + 1L
    stop("gRim must not be attempted")
  }, .package = "gRim")
  f <- fit_loglinear_graph(dat, g, engine = "loglin")
  expect_true(all(vapply(f$local_fits, function(z) z$engine == "loglin", logical(1))))
  expect_true(all(vapply(f$local_fits, function(z) z$requested_engine == "loglin", logical(1))))
  expect_true(all(vapply(f$local_fits, function(z) identical(z$attempted_engines, "loglin"), logical(1))))
  expect_true(all(vapply(f$local_fits, function(z) !isTRUE(z$fallback), logical(1))))
  expect_true(all(vapply(f$local_fits, function(z) is.na(z$fallback_reason), logical(1))))
  expect_identical(grim_calls, 0L)
})

test_that("An explicit gRim request fails instead of falling back to loglin", {
  g <- llaudit_graph("cycle4")
  dat <- llaudit_data(g)
  loglin_calls <- 0L
  local_mocked_bindings(capture.output = function(...) stop("simulated gRim crash"),
                        .package = "utils")
  local_mocked_bindings(loglin = function(...) {
    loglin_calls <<- loglin_calls + 1L
    stop("loglin must not be attempted")
  }, .package = "stats")
  e <- tryCatch(fit_loglinear_graph(dat, g, engine = "gRim"), error = identity)
  expect_s3_class(e, "atom_loglinear_error")
  expect_identical(e$engine, "gRim")
  expect_identical(e$status, "solver_failed")
  expect_identical(e$atom_id, 1L)
  expect_identical(e$attempted_engines, "gRim")
  expect_identical(names(e$engine_errors), "gRim")
  expect_match(conditionMessage(e), "simulated gRim crash", fixed = TRUE)
  expect_identical(loglin_calls, 0L)
})

test_that("A failing explicit loglin request fails without trying gRim", {
  g <- llaudit_graph("cycle4")
  dat <- llaudit_data(g)
  calls <- 0L
  local_mocked_bindings(loglin = function(...) {
    calls <<- calls + 1L
    stop("loglin down")
  }, .package = "stats")
  e <- tryCatch(fit_loglinear_graph(dat, g, engine = "loglin"), error = identity)
  expect_s3_class(e, "atom_loglinear_error")
  expect_identical(e$engine, "loglin")
  expect_identical(e$attempted_engines, "loglin")
  expect_identical(calls, 1L)
})

test_that("The removed auto engine is rejected by argument matching", {
  g <- llaudit_graph("cycle4")
  dat <- llaudit_data(g)
  expect_error(fit_loglinear_graph(dat, g, engine = "auto"), "arg")
  expect_identical(eval(formals(fit_loglinear_graph)$engine)[1L], "gRim")
})

test_that("Local validation statuses distinguish sum, dimension and margin violations", {
  g <- llaudit_graph("cycle4")
  dat <- llaudit_data(g)
  dec <- atom_decompose(g)
  ss <- margin_counts(dat, c(atoms(dec), separators(dec)))
  observed <- ss$margins[[1]]
  generators <- list(seq_along(dim(observed)))
  n <- sum(observed)
  status_of <- function(x) tryCatch({
    get(".validate_loglinear_fit", asNamespace("gAtoms"))(x, observed, generators, n, 1e-6)
    "ok"
  }, error = function(e) e$status)
  # wrong total
  expect_identical(status_of(observed * 0), "invalid_counts")
  # right total, violated generator margin
  tweaked <- array(n / length(observed), dim = dim(observed), dimnames = dimnames(observed))
  tweaked[1] <- tweaked[1] + 1
  tweaked[2] <- tweaked[2] - 1
  expect_identical(status_of(tweaked), "invalid_generator_margins")
  # wrong dimnames
  expect_identical(status_of(unname(observed)), "invalid_counts")
  # non-finite
  broken <- observed
  broken[1] <- NaN
  expect_identical(status_of(broken), "invalid_counts")
})

test_that("Explicit engine fits are deterministic across repeated identical calls", {
  g <- llaudit_graph("cycle5")
  dat <- llaudit_data(g)
  runs <- lapply(seq_len(3L), function(i) fit_loglinear_graph(dat, g, engine = "loglin"))
  engines <- lapply(runs, function(f) vapply(f$local_fits, function(z) z$engine, ""))
  expect_true(all(vapply(engines, function(e) identical(e, engines[[1L]]), logical(1))))
  expect_equal(as.numeric(logLik(runs[[1L]])), as.numeric(logLik(runs[[2L]])), tolerance = 1e-12)
  expect_equal(fitted(runs[[1L]]), fitted(runs[[2L]]), tolerance = 1e-12)
})
