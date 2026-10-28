#' Plot Atom Tree with Atom Membership (Internal)
#'
#' Core plotting logic for atom decomposition trees.
#' Internal function called by both S3 method and convenience wrapper.
#'
#' @keywords internal
#' @param atoms List whose elements are integer vectors of member vertex indices.
#' @param tree Integer matrix with two columns with dimnames c("from", "to"). 1-based indices.
#' @param labels Character vector of vertex names (optional). If NULL, uses integers.
#' @param n_vertices Total number of original graph vertices.
#' @param n_edges Total number of original graph edges.
#' @param largest_atom_size Size of the largest atom (for title/info).
#' @param vertex.size Numeric. Base node size for atom nodes. Default: 30.
#' @param vertex.label.cex Numeric. Label font size. Default: 0.8.
#' @param edge.width Numeric. Edge line width. Default: 2.
#' @param main Character. Plot title. If NULL, auto-generated.
#' @param layout Character layout mode: \code{"tree"} (default) or \code{"fr"}.
#' @param ... Additional arguments passed to \code{igraph::plot.igraph()}.
#'
#' @return NULL (invisibly). Function is called for the plot side effect.
.plot_decomp_impl <- function(atoms,
                               tree,
                               labels = NULL,
                               n_vertices = NULL,
                               n_edges = NULL,
                               largest_atom_size = NULL,
                               vertex.size = 30,
                               vertex.label.cex = 0.8,
                               edge.width = 2,
                               main = NULL,
                               layout = c("tree", "fr"),
                               ...) {

  if (!is.list(atoms) || length(atoms) == 0L) {
    stop("atoms must be a non-empty list")
  }

  layout <- match.arg(layout)

  # Normalize tree to 2-column integer matrix (possibly zero rows)
  tree <- as.matrix(tree)
  if (length(tree) == 0L) {
    tree <- matrix(integer(), ncol = 2)
  }
  if (!is.matrix(tree) || ncol(tree) != 2L) {
    stop("tree must be a matrix with 2 columns: from, to")
  }
  storage.mode(tree) <- "integer"

  # --- Extract edge list ---
  edges <- tree
  n_atoms <- length(atoms)

  if (nrow(edges) > 0L) {
    invalid <- (edges < 1L) | (edges > n_atoms)
    if (any(invalid, na.rm = TRUE)) {
      stop("tree contains atom indices outside [1, length(atoms)]")
    }
  }

  # --- Build graph with all atom vertices ---
  g <- igraph::make_empty_graph(n = n_atoms, directed = FALSE)
  if (nrow(edges) > 0L) {
    g <- igraph::add_edges(g, as.vector(t(edges)))
  }

  # --- Build node labels from atom membership ---
  # Use atoms as-is (either integer indices or character names)
  # If atoms are already character names, just join them directly
  node_labels <- vapply(seq_len(n_atoms), function(i) {
    atom_members <- atoms[[i]]
    if (is.character(atom_members)) {
      # atoms are already mapped to names
      paste(atom_members, collapse = ",")
    } else if (!is.null(labels)) {
      # atoms are integer indices, map to labels
      paste(labels[atom_members], collapse = ",")
    } else {
      # atoms are integer indices, use them directly
      paste(atom_members, collapse = ",")
    }
  }, character(1))

  # Keep node sizes in a stable visual range.
  atom_sizes <- vapply(atoms, length, integer(1))
  size_scale <- (atom_sizes - min(atom_sizes)) / max(1, diff(range(atom_sizes)))
  vertex_sizes <- vertex.size + 8 * size_scale

  # --- Auto-generate title if not provided ---
  if (is.null(main)) {
    if (n_atoms == 1L) {
      main <- "Atom Tree: Trivial decomposition (1 atom)"
    } else {
      ratio_str <- if (!is.null(largest_atom_size) && !is.null(n_vertices) && n_vertices > 0L) {
        sprintf("%.2f", largest_atom_size / n_vertices)
      } else {
        "N/A"
      }
      main <- sprintf("Atom Tree: %d atoms, max size %s (ratio: %s)",
                      n_atoms,
                      if (!is.null(largest_atom_size)) as.character(largest_atom_size) else "?",
                      ratio_str)
    }
  }

  # --- Graph layout ---
  if (layout == "tree" && igraph::ecount(g) > 0L) {
    deg <- igraph::degree(g)
    roots <- which(deg == 1L)
    if (length(roots) == 0L) {
      roots <- 1L
    }
    layout_coords <- igraph::layout_as_tree(g, root = roots[1L])
  } else {
    layout_coords <- igraph::layout_with_fr(g)
  }

  # Draw directly with igraph defaults; only enforce rectangle nodes and labels.
  igraph::plot.igraph(
    g,
    layout = layout_coords,
    vertex.shape = "rectangle",
    vertex.size = vertex_sizes,
    vertex.label = node_labels,
    vertex.label.cex = vertex.label.cex,
    edge.width = edge.width,
    margin = 0.1,
    main = main,
    ...
  )

  invisible(NULL)
}

#' Plot method for atom_decomp objects
#'
#' S3 plot method for objects returned by \code{atom_decompose()}.
#' Displays the atom tree with decomposition metadata in the title.
#'
#' @param x An object of class \code{"atom_decomp"} returned by \code{atom_decompose()}.
#' @param labels Character vector of vertex names (optional). If NULL, uses integers.
#' @param vertex.size Numeric. Base node size for atom nodes. Default: 30.
#' @param vertex.label.cex Numeric. Label font size. Default: 0.8.
#' @param edge.width Numeric. Edge line width. Default: 2.
#' @param main Character. Plot title. If NULL, auto-generated from decomposition metadata.
#' @param layout Character layout mode: \code{"tree"} (default) or \code{"fr"}.
#' @param ... Additional arguments passed to \code{igraph::plot.igraph()}.
#'
#' @return NULL (invisibly). Function is called for the plot side effect.
#'
#' @method plot atom_decomp
#' @export
plot.atom_decomp <- function(x,
                              labels = NULL,
                              vertex.size = 30,
                              vertex.label.cex = 0.8,
                              edge.width = 2,
                              main = NULL,
                              layout = c("tree", "fr"),
                              ...) {

  if (!inherits(x, "atom_decomp")) {
    stop("`x` must be an object of class 'atom_decomp'. Got class: ", paste(class(x), collapse = ", "))
  }

  if (!length(x$atoms)) {
    graphics::plot.new()
    graphics::title(main = if (is.null(main)) "Empty atom forest" else main)
    return(invisible(NULL))
  }
  dots <- list(...)
  if (is.null(dots$edge.label)) dots$edge.label <- igraph::edge_attr(atom_tree(x), "label")
  do.call(.plot_decomp_impl, c(list(atoms = x$atoms, tree = x$tree,
    labels = labels, n_vertices = igraph::vcount(x$graph),
    largest_atom_size = max(lengths(x$atoms)), vertex.size = vertex.size,
    vertex.label.cex = vertex.label.cex, edge.width = edge.width,
    main = main, layout = match.arg(layout)), dots))
  invisible(NULL)
}
