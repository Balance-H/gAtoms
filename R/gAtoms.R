#' gAtoms: structural graph decomposition
#' @useDynLib gAtoms, .registration = TRUE
"_PACKAGE"

normalize_graph <- function(graph) {
  if (inherits(graph, "atom_decomp")) graph <- graph$graph
  if (is.matrix(graph) || inherits(graph, "dgCMatrix")) {
    if (length(dim(graph)) != 2L || nrow(graph) != ncol(graph))
      stop("Adjacency matrix must be square.", call. = FALSE)
    if (is.matrix(graph) && !(is.numeric(graph) || is.logical(graph)))
      stop("Adjacency matrix must be numeric or logical.", call. = FALSE)
    values <- if (inherits(graph, "dgCMatrix")) graph@x else graph
    if (anyNA(values) || any(!is.finite(values))) stop("NA/non-finite graph structure.", call. = FALSE)
    if (any(!(values %in% c(0, 1)))) stop("Adjacency entries must be 0 or 1 (multiple/weighted edges unsupported).", call. = FALSE)
    symmetric <- if (inherits(graph, "dgCMatrix")) Matrix::isSymmetric(graph, tol = 0, check.attributes = FALSE) else isSymmetric(graph, tol = 0, check.attributes = FALSE)
    if (!isTRUE(symmetric)) stop("Adjacency matrix must be symmetric.", call. = FALSE)
    diagonal <- if (inherits(graph, "dgCMatrix")) Matrix::diag(graph) else diag(graph)
    if (any(diagonal != 0)) stop("Self-loops are unsupported.", call. = FALSE)
    rn <- rownames(graph); cn <- colnames(graph)
    if (!is.null(rn) && !is.null(cn) && !identical(rn, cn)) stop("Row and column vertex names must agree.", call. = FALSE)
    labels <- if (!is.null(rn)) rn else cn
    dimnames(graph) <- NULL
    graph <- igraph::graph_from_adjacency_matrix(graph, mode = "undirected", weighted = NULL, diag = FALSE)
    if (!is.null(labels)) graph <- igraph::set_vertex_attr(graph, "name", value = labels)
  }
  if (!igraph::is_igraph(graph)) stop("Expected an igraph, adjacency matrix or dgCMatrix.", call. = FALSE)
  if (igraph::is_directed(graph)) stop("Directed graphs are unsupported.", call. = FALSE)
  if (any(igraph::which_loop(graph))) stop("Self-loops are unsupported.", call. = FALSE)
  if (any(igraph::which_multiple(graph))) stop("Multiple edges are unsupported.", call. = FALSE)
  n <- igraph::vcount(graph)
  labels <- igraph::vertex_attr(graph, "name")
  named <- !is.null(labels)
  if (named && (length(labels) != n || anyNA(labels) || any(!nzchar(as.character(labels))) || anyDuplicated(labels)))
    stop("Vertex names must be nonempty, nonmissing and unique (no duplicated names).", call. = FALSE)
  if (named) {
    labels <- as.character(labels)
    graph <- igraph::set_vertex_attr(graph, "name", value = labels)
  } else labels <- as.character(seq_len(n))
  list(graph = graph, vertex_map = data.frame(internal_id = seq_len(n) - 1L,
       index = seq_len(n), label = labels, stringsAsFactors = FALSE), named = named)
}
.validate_graph <- function(graph) { normalize_graph(graph); invisible(TRUE) }
.external_sets <- function(sets, norm) {
  if (norm$named) lapply(sets, function(s) norm$vertex_map$label[s]) else sets
}
.prepare_graph_inputs <- function(graph) {
  e <- igraph::as_edgelist(graph, names = FALSE)
  storage.mode(e) <- "integer"
  list(edge_mat = e - 1L, directed = FALSE, n_vertices = as.integer(igraph::vcount(graph)))
}
.native_decompose <- function(graph) {
  z <- .prepare_graph_inputs(graph)
  .Call(`_gAtoms_decompose_atoms_impl`, z$edge_mat, FALSE, z$n_vertices)
}

#' Compute a minimal collapsible set using CMSA
#' @param graph An undirected simple igraph, binary adjacency matrix or dgCMatrix.
#' @param r_nodes Seed indices (one-based whole numbers), or vertex names.
#' @return One-based integer indices, named with original vertex labels when present.
#' @export
get_minimal_collapsible <- function(graph, r_nodes) {
  norm <- normalize_graph(graph); graph <- norm$graph
  if (is.character(r_nodes)) r_nodes <- match(r_nodes, norm$vertex_map$label)
  n <- igraph::vcount(graph)
  if (!is.numeric(r_nodes) || !length(r_nodes) || anyNA(r_nodes) ||
      any(!is.finite(r_nodes)) || any(r_nodes != trunc(r_nodes)) || any(r_nodes < 1 | r_nodes > n))
    stop("r_nodes must contain valid nonempty vertex names or whole-number indices.", call. = FALSE)
  z <- .prepare_graph_inputs(graph)
  out <- .Call(`_gAtoms_get_minimal_collapsible_impl`, z$edge_mat, FALSE, z$n_vertices, as.integer(unique(r_nodes)) - 1L)
  if (norm$named) names(out) <- norm$vertex_map$label[out]
  out
}

#' Expert interface to native RDA atom decomposition
#' @param graph An undirected simple igraph, binary adjacency matrix or dgCMatrix.
#' @return A named list retaining native-order `atoms`, raw RDA `separators`,
#'   and a two-column `tree` edge matrix (one-based atom IDs). `prime_ordering`
#'   is the list P1, S2, P2, ..., Sk, Pk for a connected graph; it is NULL
#'   for disconnected graphs, which have separate `component_orderings`.
#'   Vertex sets contain original names, or one-based indices for unnamed input.
#'   Fields ending in `_index` always contain one-based vertex indices.
#'   Summary metadata, graph, component membership and call are also retained.
#' @details Raw separators need not align with final tree edges. Use
#'   [atom_decompose()] for ordered separators and downstream structural work.
#' @export
decompose_atoms <- function(graph) {
  norm <- normalize_graph(graph); graph <- norm$graph
  membership <- as.integer(igraph::components(graph)$membership)
  groups <- unname(split(seq_len(igraph::vcount(graph)), membership))
  aa <- ss <- oo <- list(); edges <- matrix(integer(), ncol = 2L)
  atom_component <- integer()
  for (ci in seq_along(groups)) {
    ids <- groups[[ci]]
    z <- .native_decompose(igraph::induced_subgraph(graph, ids))
    offset <- length(aa)
    aa <- c(aa, lapply(z$atoms, function(v) ids[v]))
    ss <- c(ss, lapply(z$separators, function(v) ids[v]))
    oo[[ci]] <- lapply(z$prime_ordering, function(v) ids[v])
    if (length(z$tree)) edges <- rbind(edges, matrix(z$tree, ncol = 2L, byrow = TRUE) + offset)
    atom_component <- c(atom_component, rep.int(ci, length(z$atoms)))
  }
  colnames(edges) <- c("from", "to")
  out <- list(atoms = .external_sets(aa, norm), separators = .external_sets(ss, norm), tree = edges,
    prime_ordering = if (length(oo) == 1L) .external_sets(oo[[1L]], norm) else if (!length(oo)) list() else NULL,
    component_orderings = lapply(oo, .external_sets, norm = norm), atoms_index = aa, separators_index = ss,
    component_orderings_index = oo, graph = graph, vertex_map = norm$vertex_map,
    components = list(vertices = membership, atoms = atom_component), n_vertices = igraph::vcount(graph),
    n_edges = igraph::ecount(graph), n_atoms = length(aa), largest_atom_size = max(c(0L, lengths(aa))),
    median_atom_size = if (length(aa)) stats::median(lengths(aa)) else NA_real_,
    largest_separator_size = max(c(0L, lengths(ss))), call = match.call())
  out$largest_atom_ratio <- if (out$n_vertices) out$largest_atom_size / out$n_vertices else NA_real_
  out
}

#' Convert supported graph representations
#' @param x An atom_decomp object or any supported graph input.
#' @return `as_igraph()` returns an igraph; `as_dgCMatrix()` returns a sparse
#'   binary adjacency matrix. Original labels are preserved. Unnamed vertices
#'   retain deterministic one-based positional indices.
#' @export
as_igraph <- function(x) normalize_graph(x)$graph
#' @rdname as_igraph
#' @export
as_dgCMatrix <- function(x) {
  g <- as_igraph(x)
  methods::as(igraph::as_adjacency_matrix(g, sparse = TRUE), "dgCMatrix")
}
