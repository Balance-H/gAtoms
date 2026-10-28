# Audit-only tests: local generator construction for graphical log-linear atoms.
# Reference cliques are computed by an independent brute-force enumeration and,
# where available, cross-checked against gRbase; the production helper
# .loglinear_generators is never used as its own reference.

llaudit_graph <- function(kind) {
  g <- switch(kind,
    path4 = igraph::make_ring(4, circular = FALSE),
    path5 = igraph::make_ring(5, circular = FALSE),
    cycle4 = igraph::make_ring(4),
    cycle5 = igraph::make_ring(5),
    complete3 = igraph::make_full_graph(3),
    complete4 = igraph::make_full_graph(4),
    mixed = igraph::add_edges(igraph::add_vertices(igraph::make_ring(4), 2),
                              c(3, 5, 4, 5, 5, 6)),
    twocliques = igraph::graph_from_edgelist(
      rbind(c(1, 2), c(2, 3), c(1, 3), c(3, 4), c(4, 5), c(3, 5)), directed = FALSE),
    disco = igraph::disjoint_union(igraph::make_ring(3, circular = FALSE),
                                   igraph::make_empty_graph(1, directed = FALSE)),
    disco2 = igraph::disjoint_union(igraph::make_ring(4),
                                    igraph::make_full_graph(2)))
  igraph::set_vertex_attr(g, "name", value = LETTERS[seq_len(igraph::vcount(g))])
}

llaudit_data <- function(g) {
  grid <- do.call(expand.grid,
    stats::setNames(rep(list(c("lo", "hi")), igraph::vcount(g)), igraph::V(g)$name))
  reps <- rep(3L, nrow(grid))
  reps[seq.int(4L, 7L)] <- 4L
  grid[rep(seq_len(nrow(grid)), reps), , drop = FALSE]
}

# Independent brute-force maximal-clique enumeration for a small igraph.
llaudit_max_cliques_bruteforce <- function(g) {
  n <- igraph::vcount(g)
  adj <- igraph::as_adjacency_matrix(g, sparse = FALSE) != 0
  is_clique <- function(s) {
    if (length(s) < 2L) return(TRUE)
    sub <- adj[s, s, drop = FALSE]
    all(sub[upper.tri(sub)])
  }
  cliques <- list()
  for (m in seq_len(n)) {
    cmb <- utils::combn(n, m, simplify = FALSE)
    for (s in cmb) if (is_clique(s)) cliques[[length(cliques) + 1L]] <- sort(s)
  }
  maximal <- vapply(cliques, function(s) !any(vapply(cliques, function(t)
    length(t) > length(s) && all(s %in% t), logical(1))), logical(1))
  keys <- vapply(cliques[maximal], function(s) paste(s, collapse = ","), "")
  cliques[maximal][order(keys, method = "radix")]
}

llaudit_key <- function(s) paste(sort(as.integer(s)), collapse = ",")

llaudit_gen_sets <- function(z) {
  # generators stored in local induced-subgraph integer indices
  lapply(z$generators, function(v) sort(match(v, z$vertices)))
}

test_that("Generators equal the maximal cliques of the induced atom graph G[P_i]", {
  for (kind in c("complete3", "complete4", "path4", "path5", "cycle4", "cycle5",
                 "mixed", "twocliques", "disco", "disco2")) {
    g <- llaudit_graph(kind)
    dec <- atom_decompose(g)
    for (i in seq_along(dec$atoms_index)) {
      v <- sort(dec$atoms_index[[i]])
      sub <- igraph::induced_subgraph(g, v)
      reference <- llaudit_max_cliques_bruteforce(sub)
      if (kind %in% c("complete4", "path5", "cycle5")) {
        # cross-check the brute-force oracle itself against igraph
        expect_setequal(lapply(unique(reference), llaudit_key),
                        lapply(unique(lapply(igraph::max_cliques(sub), sort)), llaudit_key))
      }
      if (dec$complete[i]) {
        expect_length(reference, 1L)
        expect_identical(reference[[1L]], seq_along(v))
      }
      got <- get(".loglinear_generators", asNamespace("gAtoms"))(sub)
      expect_setequal(lapply(got, llaudit_key), lapply(reference, llaudit_key))
    }
  }
})

test_that("Generators returned by a fitted object are the induced maximal cliques in local indices", {
  for (kind in c("path4", "cycle4", "mixed", "twocliques", "disco")) {
    g <- llaudit_graph(kind)
    dat <- llaudit_data(g)
    dec <- atom_decompose(g)
    f <- fit_loglinear_graph(dat, g, engine = "loglin")
    for (i in seq_along(f$local_fits)) {
      z <- f$local_fits[[i]]
      sub <- igraph::induced_subgraph(g, sort(dec$atoms_index[[i]]))
      reference <- llaudit_max_cliques_bruteforce(sub)
      expect_identical(z$vertices, igraph::V(sub)$name)
      expect_setequal(lapply(z$generators, function(x) paste(sort(x), collapse = ",")),
                      lapply(reference, function(s) paste(igraph::V(sub)$name[s], collapse = ",")))
      if (dec$complete[i]) {
        expect_identical(z$engine, "closed_form")
        expect_identical(z$generators, list(z$vertices))
      }
    }
  }
})

test_that("Generator indices are valid, clique, maximal, duplicate-free and ordered", {
  for (kind in c("path4", "cycle4", "cycle5", "mixed", "twocliques", "disco2")) {
    g <- llaudit_graph(kind)
    dec <- atom_decompose(g)
    for (i in seq_along(dec$atoms_index)) {
      sub <- igraph::induced_subgraph(g, sort(dec$atoms_index[[i]]))
      got <- get(".loglinear_generators", asNamespace("gAtoms"))(sub)
      q <- igraph::vcount(sub)
      expect_true(all(vapply(got, function(s) is.integer(s) && length(s) > 0L &&
        !anyNA(s) && all(s >= 1L & s <= q) && !anyDuplicated(s), logical(1))))
      # sorted ascending
      expect_true(all(vapply(got, function(s) identical(s, sort(s)), logical(1))))
      # every generator is a clique of G[P_i]
      expect_true(all(vapply(got, function(s) gAtoms:::.is_complete(sub, s), logical(1))))
      # maximal: no generator is a proper subset of another
      expect_false(any(vapply(seq_along(got), function(a) any(vapply(got[-a],
        function(s) length(s) < length(got[[a]]) && all(s %in% got[[a]]), logical(1))), logical(1))))
      # unique
      expect_false(anyDuplicated(vapply(got, llaudit_key, "")) > 0L)
      # deterministic radix order of the keys
      keys <- vapply(got, llaudit_key, "")
      expect_identical(keys, sort(keys, method = "radix"))
      # no missing maximal clique: every vertex-time-edge pair is covered
      reference <- llaudit_max_cliques_bruteforce(sub)
      expect_setequal(lapply(got, llaudit_key), lapply(reference, llaudit_key))
    }
  }
})

test_that("Nameless vertices still yield induced maximal cliques", {
  g <- igraph::add_edges(igraph::make_empty_graph(4, directed = FALSE),
                         c(1, 2, 2, 3, 3, 1, 3, 4))
  # normalize_graph assigns positional labels for unnamed input
  sub <- igraph::induced_subgraph(g, 1:4)
  got <- get(".loglinear_generators", asNamespace("gAtoms"))(sub)
  reference <- llaudit_max_cliques_bruteforce(sub)
  expect_setequal(lapply(got, llaudit_key), lapply(reference, llaudit_key))
  expect_identical(got, list(c(1L, 2L, 3L), c(3L, 4L)))
})

test_that("gRbase independently agrees with the induced maximal cliques", {
  skip_if_not_installed("gRbase")
  for (kind in c("cycle4", "cycle5", "mixed", "twocliques")) {
    g <- llaudit_graph(kind)
    dec <- atom_decompose(g)
    for (i in seq_along(dec$atoms_index)) {
      sub <- igraph::induced_subgraph(g, sort(dec$atoms_index[[i]]))
      got <- get(".loglinear_generators", asNamespace("gAtoms"))(sub)
      amat <- igraph::as_adjacency_matrix(sub, sparse = FALSE)
      amat <- amat[order(igraph::V(sub)$name), order(igraph::V(sub)$name), drop = FALSE]
      ref <- gRbase::max_cliqueMAT(amat)
      refkeys <- lapply(ref, function(s) llaudit_key(match(s, colnames(amat))))
      expect_setequal(lapply(got, llaudit_key), refkeys)
    }
  }
})

test_that("Both local engines receive exactly the package generators", {
  g <- llaudit_graph("cycle4")
  dat <- llaudit_data(g)
  dec <- atom_decompose(g)
  expected <- get(".loglinear_generators", asNamespace("gAtoms"))(
    igraph::induced_subgraph(g, sort(dec$atoms_index[[1]])))
  seen <- NULL
  real_loglin <- get("loglin", envir = asNamespace("stats"))
  local_mocked_bindings(loglin = function(table, margin, ...) {
    seen <<- margin
    real_loglin(table = table, margin = margin, ...)
  }, .package = "stats")
  f <- fit_loglinear_graph(dat, g, engine = "loglin")
  expect_identical(seen, expected)
  expect_identical(lapply(f$local_fits[[1]]$generators,
                          function(v) sort(match(v, f$local_fits[[1]]$vertices))), expected)
})

test_that("gRim adapter passes one margin per generator, in order", {
  skip_if_not_installed("gRim")
  g <- llaudit_graph("cycle4")
  dat <- llaudit_data(g)
  f <- fit_loglinear_graph(dat, g, engine = "gRim")
  z <- f$local_fits[[1]]
  expect_identical(z$engine, "gRim")
  expect_identical(z$method, "effloglin")
  expect_identical(dim(z$fitted_counts), dim(z$observed))
  expect_identical(dimnames(z$fitted_counts), dimnames(z$observed))
  # the fitted counts satisfy exactly the generator margins, i.e. the margins the
  # adapter handed to gRim were these generators
  for (gen in z$generators) {
    s <- match(gen, z$vertices)
    expect_equal(as.vector(margin.table(z$fitted_counts, s)),
                 as.vector(margin.table(z$observed, s)), tolerance = 1e-6)
  }
})
