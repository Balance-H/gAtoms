.forest_graph <- function(k, edges) {
  g <- igraph::make_empty_graph(k, directed = FALSE)
  if (length(edges)) g <- igraph::add_edges(g, as.vector(t(edges)))
  g
}
.is_complete <- function(g, v) {
  length(v) < 2L || igraph::ecount(igraph::induced_subgraph(g, v)) == length(v) * (length(v) - 1) / 2
}
.set_key <- function(v) paste(sort(v), collapse = ",")

#' Construct an ordered structural atom decomposition
#' @param graph An undirected simple igraph, binary adjacency matrix or dgCMatrix.
#' @return An S3 `atom_decomp` list with `graph`, `vertex_map`, ordered `atoms`,
#'   ordered `separators` (one per nonroot atom), `separator_atoms` (their atom
#'   IDs), `residuals`, `parent` (NA at roots), `children`, `tree` (edge matrix),
#'   `trees` (component igraph trees), `components` (vertex/atom memberships),
#'   `complete`, `native_atom_ids`, `prime_ordering`, `component_orderings`,
#'   `diagnostics`, `call`, and `version`. Vertex-set `_index` fields use R's
#'   one-based positions; `vertex_map$internal_id` documents C's zero-based IDs.
#'   Raw native separators are retained in `diagnostics$raw_separators_index`.
#' @details Within each component, the ordering is P1, S2, P2, ..., Sk, Pk.
#'   Si is Pi intersected with its earlier tree parent and with the union of
#'   previous atoms. Residuals are P1 at the root and Pi minus Si elsewhere.
#'   All atom indices refer to the final ordered atoms. Each graph component
#'   has its own tree and ordering; no artificial cross-component edges or
#'   empty separators are introduced. Empty graphs have no atoms or orderings.
#' @examples
#' dec <- atom_decompose(igraph::make_ring(6))
#' atoms(dec)
#' prime_ordering(dec)
#' validate(dec)
#' @export
atom_decompose <- function(graph) {
  start <- proc.time()[[3L]]
  raw <- decompose_atoms(graph); norm <- normalize_graph(raw$graph)
  oo <- raw$component_orderings_index
  aa <- unlist(lapply(oo, function(o) o[seq.int(1L, length(o), by = 2L)]), recursive = FALSE)
  if (is.null(aa)) aa <- list()
  native <- match(vapply(aa, .set_key, ""), vapply(raw$atoms_index, .set_key, ""))
  inverse <- match(seq_along(raw$atoms_index), native)
  tree <- raw$tree
  if (length(tree)) tree[] <- inverse[tree]
  k <- length(aa); parent <- rep.int(NA_integer_, k)
  for (i in seq_len(nrow(tree))) {
    child <- max(tree[i, ]); par <- min(tree[i, ])
    if (!is.na(parent[child])) stop("Invalid native tree: multiple earlier neighbours.", call. = FALSE)
    parent[child] <- par
  }
  nonroots <- which(!is.na(parent))
  ss <- lapply(nonroots, function(i) intersect(aa[[i]], aa[[parent[i]]]))
  rr <- lapply(seq_len(k), function(i) if (is.na(parent[i])) aa[[i]] else setdiff(aa[[i]], aa[[parent[i]]]))
  ac <- raw$components$atoms[native]
  trees <- lapply(seq_along(oo), function(ci) {
    ids <- which(ac == ci); e <- tree[tree[, 1L] %in% ids, , drop = FALSE]
    if (length(e)) e[] <- match(e, ids)
    igraph::set_vertex_attr(.forest_graph(length(ids), e), "atom_id", value = ids)
  })
  out <- list(graph = raw$graph, vertex_map = raw$vertex_map,
    atoms = .external_sets(aa, norm), separators = .external_sets(ss, norm),
    prime_ordering = raw$prime_ordering, component_orderings = raw$component_orderings,
    residuals = .external_sets(rr, norm), parent = parent,
    children = lapply(seq_len(k), function(i) which(parent == i)), tree = tree, trees = trees,
    components = list(vertices = raw$components$vertices, atoms = ac),
    complete = vapply(aa, function(v) .is_complete(raw$graph, v), logical(1)),
    native_atom_ids = native, separator_atoms = nonroots,
    atoms_index = aa, separators_index = ss, residuals_index = rr,
    component_orderings_index = oo, diagnostics = list(raw_separators_index = raw$separators_index,
      elapsed = unname(proc.time()[[3L]] - start)), call = match.call(), version = "0.1.0")
  class(out) <- "atom_decomp"
  out
}

#' Access structural decomposition fields
#' @param x An atom_decomp object.
#' @param ... Additional method arguments.
#' @return `atoms`, `separators`, and `residuals` return vertex-set lists.
#'   `atom_tree` returns the igraph atom forest with ordered atom IDs and
#'   separator edge labels. `prime_ordering` returns an alternating list for
#'   connected input and a list of component orderings for disconnected input.
#' @export
atoms <- function(x, ...) UseMethod("atoms")
#' @rdname atoms
#' @export
atoms.atom_decomp <- function(x, ...) x$atoms
#' @rdname atoms
#' @export
separators <- function(x, ...) UseMethod("separators")
#' @rdname atoms
#' @export
separators.atom_decomp <- function(x, ...) x$separators
#' @rdname atoms
#' @param object An atom_decomp object (for the standard stats residuals generic).
#' @export
residuals.atom_decomp <- function(object, ...) object$residuals
#' @rdname atoms
#' @export
prime_ordering <- function(x, ...) UseMethod("prime_ordering")
#' @rdname atoms
#' @export
prime_ordering.atom_decomp <- function(x, ...) {
  if (length(x$component_orderings) <= 1L) x$prime_ordering else x$component_orderings
}
#' @rdname atoms
#' @export
atom_tree <- function(x, ...) UseMethod("atom_tree")
#' @rdname atoms
#' @export
atom_tree.atom_decomp <- function(x, ...) {
  g <- .forest_graph(length(x$atoms), x$tree)
  g <- igraph::set_vertex_attr(g, "atom_id", value = seq_along(x$atoms))
  labels <- vapply(seq_len(nrow(x$tree)), function(i) paste(intersect(x$atoms[[x$tree[i, 1L]]], x$atoms[[x$tree[i, 2L]]]), collapse = ","), "")
  igraph::set_edge_attr(g, "label", value = labels)
}

#' Structural decomposition cost diagnostics
#' @param x An atom_decomp object.
#' @return A named list of task-independent graph and atom size diagnostics.
#' @export
decomp_cost <- function(x) {
  p <- igraph::vcount(x$graph); k <- length(x$atoms)
  list(p = p, m = igraph::ecount(x$graph), n_components = length(x$component_orderings),
    n_atoms = k, max_atom_size = max(c(0L, lengths(x$atoms))),
    max_separator_size = max(c(0L, lengths(x$separators))), sum_atom_sizes = sum(lengths(x$atoms)),
    largest_prime_ratio = if (p) max(lengths(x$atoms)) / p else NA_real_,
    complete_atom_count = sum(x$complete), complete_atom_proportion = if (k) mean(x$complete) else NA_real_)
}
#' Print an atom decomposition
#' @param x An atom_decomp object.
#' @param ... Unused.
#' @return x invisibly.
#' @method print atom_decomp
#' @export
print.atom_decomp <- function(x, ...) {
  d <- decomp_cost(x)
  cat(sprintf("Atom decomposition: %d vertices, %d edges, %d components, %d atoms\n", d$p, d$m, d$n_components, d$n_atoms))
  cat(sprintf("Maximum atom: %d; maximum separator: %d; largest-prime ratio: %.3f\n", d$max_atom_size, d$max_separator_size, d$largest_prime_ratio))
  if (!is.null(x$diagnostics$elapsed)) cat(sprintf("Decomposition time: %.3f seconds\n", x$diagnostics$elapsed))
  invisible(x)
}
#' Summarize an atom decomposition
#' @param object An atom_decomp object.
#' @param ... Unused.
#' @return A structured list with size distributions, component sizes and cost diagnostics.
#' @method summary atom_decomp
#' @export
summary.atom_decomp <- function(object, ...) {
  c(decomp_cost(object), list(atom_size_distribution = table(lengths(object$atoms)),
    separator_size_distribution = table(lengths(object$separators)),
    incomplete_atom_count = sum(!object$complete), component_vertex_sizes = table(object$components$vertices),
    component_atom_sizes = table(object$components$atoms)))
}
