test_that("300 seeded small and medium graphs agree with optional gRbase mpd", {
  skip_on_cran()
  skip_if_not_installed("gRbase")
  set.seed(20260910)
  canonical <- function(sets) sort(unname(vapply(sets, function(v) paste(sort(as.integer(v)), collapse = ","), "")))
  for (i in 1:300) {
    n <- if (i <= 240) sample(2:20, 1) else sample(25:60, 1)
    g <- igraph::sample_gnp(n, sample(c(.03, .1, .25, .5, .8), 1))
    d <- atom_decompose(g)
    expect_identical(validate(d), TRUE, info = paste("case", i))
    # Compare each connected component; labels preserve its original indices.
    groups <- split(seq_len(n), igraph::components(g)$membership)
    oracle <- list()
    for (ids in groups) {
      if (length(ids) == 1L) { oracle <- c(oracle, list(ids)); next }
      sub <- igraph::induced_subgraph(g, ids)
      sub <- igraph::set_vertex_attr(sub, "name", value = as.character(ids))
      adj <- igraph::as_adjacency_matrix(sub, sparse = TRUE)
      # gRbase 2.0.3 mpd.default uses triangulate(), despite the generic's
      # minimal_triang() default. A nonminimal triangulation can merge atoms
      # (regression: seeded case 248). Specify the minimal triangulation.
      ref <- gRbase::mpd(adj, tobject = gRbase::minimal_triang(adj))
      oracle <- c(oracle, ref$cliques)
    }
    expect_identical(canonical(d$atoms_index), canonical(oracle), info = paste("case", i))
  }
})
