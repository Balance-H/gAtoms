test_that("get_minimal_collapsible returns 1-based integer vector", {
  g <- igraph::make_ring(6)
  out <- get_minimal_collapsible(g, c(1L, 3L))

  expect_type(out, "integer")
  expect_true(length(out) >= 1L)
  expect_true(all(out >= 1L & out <= igraph::vcount(g)))
})

test_that("decompose_atoms returns named list", {
  g <- igraph::make_ring(6)
  out <- decompose_atoms(g)

  expect_true(all(c("atoms", "separators", "tree", "prime_ordering") %in% names(out)))
  expect_type(out$atoms, "list")
  expect_type(out$separators, "list")
  expect_true(is.matrix(out$tree))
})

test_that("index conversion is stable for complete graph seeds", {
  g <- igraph::make_full_graph(5)

  # Full graph is already clique-separable; seed set should stay unchanged.
  out <- get_minimal_collapsible(g, c(2L, 3L))

  expect_identical(out, c(2L, 3L))
})

