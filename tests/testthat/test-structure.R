canon <- function(sets) sort(vapply(sets, function(v) paste(sort(as.character(v)), collapse = ","), ""))
native_order <- function(atoms, edges) {
  .Call(gAtoms:::`_gAtoms_build_prime_ordering_impl`, lapply(atoms, function(v) as.integer(v) - 1L),
        matrix(as.integer(edges) - 1L, ncol = 2L))
}
check_dec <- function(g) {
  d <- atom_decompose(g)
  expect_identical(validate(d, expensive = TRUE), TRUE)
  expect_identical(class(d), "atom_decomp")
  d
}
test_that("leaf peeling has exact deterministic chain and star results", {
  expect_identical(native_order(list(3:1), matrix(integer(), 0, 2)), list(3:1))
  a <- list(c(3L, 1L), c(2L, 3L))
  expect_identical(native_order(a, matrix(c(1, 2), 1)), list(a[[2]], 3L, a[[1]]))
  a <- list(c(5L, 1L, 3L), c(3L, 2L, 1L))
  expect_identical(native_order(a, matrix(c(1, 2), 1)), list(a[[2]], c(1L, 3L), a[[1]]))
  a <- list(c(2L, 1L), c(3L, 2L), c(4L, 3L), c(5L, 4L))
  e <- cbind(1:3, 2:4)
  expect_identical(native_order(a, e), list(a[[4]], 4L, a[[3]], 3L, a[[2]], 2L, a[[1]]))
  a <- list(c(9L, 1L, 2L), c(3L, 9L), c(9L, 4L), c(5L, 9L))
  e <- cbind(c(1, 1, 1), 2:4)
  expected <- list(a[[4]], 9L, a[[1]], 9L, a[[3]], 9L, a[[2]])
  expect_identical(native_order(a, e), expected)
  for (i in 1:20) expect_identical(native_order(a, e), expected)
  # Independent deep copies survive mutation of either returned atom list.
  z <- decompose_atoms(igraph::make_ring(4)); old <- z$prime_ordering
  z$atoms[[1]][1] <- 999L
  expect_identical(z$prime_ordering, old)
})
test_that("invalid native atom trees and native graph inputs fail safely", {
  a <- list(1:2, 2:3, 3:4)
  expect_error(native_order(a, matrix(c(1, 2), 1)), "atom tree")
  expect_error(native_order(a, rbind(c(1, 2), c(1, 2))), "atom tree")
  expect_error(native_order(a, rbind(c(1, 2), c(2, 3), c(3, 1))), "atom tree")
  expect_identical(native_order(list(), matrix(integer(), 0, 2)), list())
  f <- gAtoms:::`_gAtoms_decompose_atoms_impl`
  expect_error(.Call(f, matrix(c(0L, 1L), 1), TRUE, 2L), "directed")
  expect_error(.Call(f, matrix(c(-1L, 1L), 1), FALSE, 2L), "endpoints")
  expect_error(.Call(f, matrix(c(0L, 0L), 1), FALSE, 1L), "simple")
  expect_error(.Call(f, matrix(integer(), 0, 2), FALSE, 2L), "connected")
  expect_error(.Call(f, matrix(integer(), 0, 2), FALSE, 1.5), "whole")
  expect_error(.Call(f, matrix(0.5, 1, 2), FALSE, 2L), "whole")
  # Native error handlers must have been restored: subsequent calls still work.
  check_dec(igraph::make_ring(5))
})
test_that("one atom, two atoms, chains and branching decompositions validate", {
  expect_length(prime_ordering(check_dec(igraph::make_full_graph(5))), 1L)
  d <- check_dec(igraph::make_graph(c(1, 2, 2, 3), directed = FALSE))
  expect_length(prime_ordering(d), 3L)
  expect_setequal(d$separators[[1]], intersect(d$atoms[[1]], d$atoms[[2]]))
  check_dec(igraph::make_ring(7, circular = FALSE))
  check_dec(igraph::make_star(7, mode = "undirected"))
  expect_identical(check_dec(igraph::make_ring(7))$complete, FALSE)
})
test_that("supported inputs preserve labels and sparse round trips", {
  g <- igraph::set_vertex_attr(igraph::make_ring(6), "name", value = c("z", "a", "q", "2", "10", "v"))
  sparse <- as_dgCMatrix(g); expect_s4_class(sparse, "dgCMatrix")
  for (input in list(g, as.matrix(sparse), sparse)) {
    d <- check_dec(input)
    expect_identical(d$vertex_map$label, igraph::V(g)$name)
    expect_setequal(unlist(atoms(d)), igraph::V(g)$name)
    expect_identical(igraph::V(as_igraph(as_dgCMatrix(d)))$name, igraph::V(g)$name)
    expect_identical(d$vertex_map$internal_id, 0:5)
    expect_true(all(unlist(d$atoms_index) >= 1L))
  }
  d <- check_dec(igraph::make_ring(5))
  expect_identical(d$vertex_map$label, as.character(1:5))
  expect_type(d$atoms[[1]], "integer")
  expect_null(igraph::vertex_attr(as_igraph(as_dgCMatrix(d)), "name"))
  expect_error(get_minimal_collapsible(g, 1.5), "whole")
  expect_identical(unname(get_minimal_collapsible(g, c("z", "q"))), 1:6)
})
test_that("invalid structures are rejected without silent repair", {
  expect_error(atom_decompose(igraph::make_ring(4, directed = TRUE)), "Directed")
  expect_error(atom_decompose(igraph::make_graph(c(1, 1), directed = FALSE)), "Self-loop")
  expect_error(atom_decompose(igraph::make_graph(c(1, 2, 1, 2), directed = FALSE)), "Multiple")
  g <- igraph::set_vertex_attr(igraph::make_ring(4), "name", value = c("a", "a", "b", "c"))
  expect_error(atom_decompose(g), "duplicated")
  expect_error(atom_decompose(matrix(0, 2, 3)), "square")
  expect_error(atom_decompose(matrix(c(0, 1, 0, 0), 2)), "symmetric")
  expect_error(atom_decompose(matrix(NA_real_, 2, 2)), "NA")
  expect_error(atom_decompose(diag(2)), "Self-loop")
  expect_error(atom_decompose(matrix(2, 2, 2)), "0 or 1")
  m <- matrix(0, 2, 2, dimnames = list(c("a", "b"), c("b", "a")))
  expect_error(atom_decompose(m), "agree")
  m <- Matrix::sparseMatrix(i = 1L, j = 2L, x = NA_real_, dims = c(2, 2))
  expect_error(atom_decompose(m), "NA")
})
test_that("empty, disconnected and isolated graphs have separate orderings", {
  d <- check_dec(igraph::make_empty_graph(0, directed = FALSE))
  expect_length(atoms(d), 0L); expect_length(prime_ordering(d), 0L)
  expect_identical(dim(as_dgCMatrix(d)), c(0L, 0L))
  expect_length(prime_ordering(check_dec(igraph::make_empty_graph(1, directed = FALSE))), 1L)
  g <- igraph::disjoint_union(igraph::make_ring(5), igraph::make_ring(4, circular = FALSE), igraph::make_empty_graph(2, directed = FALSE))
  d <- check_dec(g)
  expect_length(d$trees, 4L); expect_length(d$component_orderings, 4L)
  expect_equal(nrow(d$tree), length(d$atoms) - 4L)
  expect_true(all(lengths(d$separators) > 0L))
  expect_null(d$prime_ordering)
  expect_equal(lengths(d$component_orderings), c(1, 5, 1, 1))
})
test_that("S3 accessors, diagnostics and plots work", {
  d <- check_dec(igraph::make_star(6, mode = "undirected"))
  expect_output(print(d), "6 vertices")
  expect_equal(summary(d)$n_atoms, length(atoms(d)))
  expect_identical(stats::residuals(d), d$residuals)
  expect_identical(separators(d), d$separators)
  expect_equal(igraph::vcount(atom_tree(d)), length(d$atoms))
  expect_equal(decomp_cost(d)$sum_atom_sizes, sum(lengths(d$atoms)))
  path <- tempfile(fileext = ".pdf"); grDevices::pdf(path)
  on.exit({grDevices::dev.off(); unlink(path)}, add = TRUE)
  expect_no_error(plot(d)); expect_no_error(plot(atom_decompose(igraph::make_empty_graph(0, directed = FALSE))))
  expect_no_error(plot(atom_decompose(igraph::disjoint_union(igraph::make_ring(4), igraph::make_ring(3)))))
})
test_that("validator detects corrupted metadata and mathematical invariants", {
  d <- atom_decompose(igraph::make_ring(5, circular = FALSE))
  for (field in c("parent", "children", "complete", "vertex_map", "component_orderings", "residuals_index", "separators_index", "native_atom_ids")) {
    bad <- d; bad[[field]] <- NULL
    result <- validate(bad)
    expect_s3_class(result, "atom_decomp_validation")
    expect_true(length(result$errors) > 0L)
  }
  bad <- d; bad$atoms_index[[1]] <- integer(); expect_type(validate(bad)$errors, "character")
  bad <- d; bad$tree[1, ] <- c(1L, 1L); expect_type(validate(bad)$errors, "character")
  bad <- d; bad$atoms_index[[length(bad$atoms_index)]] <- unique(c(bad$atoms_index[[length(bad$atoms_index)]], d$atoms_index[[1]][1]))
  expect_true(any(grepl("running intersection", validate(bad)$errors)))
})
test_that("random small graphs satisfy invariants deterministically", {
  set.seed(4102)
  for (i in 1:80) {
    g <- igraph::sample_gnp(sample(2:18, 1), runif(1, .05, .85))
    d <- check_dec(g)
    expect_identical(d$component_orderings, atom_decompose(g)$component_orderings)
  }
})

test_that("complete flags agree with an independent pairwise clique test", {
  g <- igraph::disjoint_union(igraph::make_full_graph(4), igraph::make_ring(5),
                              igraph::make_star(5, mode = "undirected"))
  d <- atom_decompose(g)
  pairwise <- vapply(d$atoms_index, function(v) {
    if (length(v) < 2L) return(TRUE)
    pairs <- utils::combn(v, 2L)
    all(vapply(seq_len(ncol(pairs)), function(i) igraph::are_adjacent(g, pairs[1L, i], pairs[2L, i]), logical(1)))
  }, logical(1))
  expect_identical(d$complete, pairwise)
  expect_true(any(d$complete) && any(!d$complete))
})
