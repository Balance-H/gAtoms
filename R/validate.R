#' Validate a structural decomposition
#' @param x An atom_decomp object.
#' @param expensive Also check graph edge coverage and atom maximality.
#' @param ... Additional method arguments.
#' @return TRUE on success, otherwise an `atom_decomp_validation` list with
#'   `valid = FALSE` and a character vector `errors` describing failed invariants.
#' @details Checks coverage, valid unique vertex sets, component membership,
#'   mapping, tree topology, edge intersections, clique separators, alternating
#'   perfect orderings, previous-union intersections, parent/children, residuals,
#'   complete flags and the running-intersection property. No reference ordering
#'   from another package is required. Default checks use sparse graph storage.
#' @export
validate <- function(x, ...) UseMethod("validate")
#' @rdname validate
#' @export
validate.atom_decomp <- function(x, expensive = FALSE, ...) {
  errors <- character()
  check <- function(ok, text) if (!isTRUE(ok)) errors <<- c(errors, text)
  tryCatch({
    required <- c("graph", "vertex_map", "atoms", "separators", "residuals", "parent", "children", "tree", "trees",
      "components", "complete", "native_atom_ids", "separator_atoms", "atoms_index", "separators_index",
      "residuals_index", "component_orderings", "component_orderings_index", "prime_ordering")
    missing <- setdiff(required, names(x))
    if (length(missing)) stop("Missing fields: ", paste(missing, collapse = ", "))
    norm <- normalize_graph(x$graph); g <- norm$graph; p <- igraph::vcount(g)
    aa <- x$atoms_index; k <- length(aa)
    check(is.list(aa), "atoms_index must be a list")
    valid_set <- function(v, empty = FALSE) is.numeric(v) && (empty || length(v) > 0L) &&
      !anyNA(v) && all(is.finite(v)) && all(v == trunc(v)) && all(v >= 1L & v <= p) && !anyDuplicated(v)
    good <- is.list(aa) && all(vapply(aa, valid_set, logical(1)))
    check(good, "invalid or empty atom vertex set")
    if (!good) stop("Cannot inspect malformed atoms")
    check(!anyDuplicated(vapply(aa, .set_key, "")), "duplicate atoms")
    check(setequal(unlist(aa), seq_len(p)), "atoms do not cover graph vertices")
    check(identical(x$vertex_map, norm$vertex_map), "internal/external vertex mapping mismatch")
    check(identical(x$atoms, .external_sets(aa, norm)), "external atom labels mismatch")
    check(identical(x$separators, .external_sets(x$separators_index, norm)), "external separator labels mismatch")
    check(identical(x$residuals, .external_sets(x$residuals_index, norm)), "external residual labels mismatch")
    check(setequal(x$native_atom_ids, seq_len(k)) && length(x$native_atom_ids) == k, "native atom ID mapping is not a permutation")
    vm <- as.integer(igraph::components(g)$membership); nc <- length(unique(vm))
    check(identical(x$components$vertices, vm), "vertex component assignments mismatch")
    ac <- vapply(aa, function(a) vm[a[1L]], integer(1))
    check(all(vapply(seq_along(aa), function(i) all(vm[aa[[i]]] == ac[i]), logical(1))), "atom crosses graph components")
    check(identical(x$components$atoms, ac), "atom component assignments mismatch")
    e <- x$tree
    good_tree <- is.matrix(e) && ncol(e) == 2L && !anyNA(e) && all(e == trunc(e)) && all(e >= 1L & e <= k)
    check(good_tree, "invalid atom tree edge matrix")
    if (!good_tree) stop("Cannot inspect malformed tree")
    forest <- .forest_graph(k, e)
    check(igraph::vcount(forest) == k, "tree vertex count mismatch")
    check(igraph::is_simple(forest), "atom forest has loops or multiple edges")
    check(nrow(e) == k - nc, "forest edge count mismatch")
    check(all(ac[e[, 1L]] == ac[e[, 2L]]), "cross-component atom tree edge")
    check(length(x$trees) == nc, "component tree count mismatch")
    for (ci in seq_len(nc)) {
      ids <- which(ac == ci); sub <- igraph::induced_subgraph(forest, ids)
      check(igraph::is_connected(sub) && igraph::ecount(sub) == length(ids) - 1L, paste("component", ci, "tree is not connected/acyclic"))
      stored <- x$trees[[ci]]
      check(igraph::vcount(stored) == length(ids) && identical(as.integer(igraph::vertex_attr(stored, "atom_id")), ids), paste("component", ci, "tree vertex IDs mismatch"))
      key_edges <- function(t) sort(apply(igraph::as_edgelist(t, names = FALSE), 1L, function(z) paste(sort(z), collapse = ",")))
      check(identical(key_edges(stored), key_edges(sub)), paste("component", ci, "stored tree edges mismatch"))
    }
    check(length(x$parent) == k && length(x$children) == k, "parent/children lengths mismatch")
    expected_parent <- rep.int(NA_integer_, k)
    for (i in seq_len(nrow(e))) {
      child <- max(e[i, ]); par <- min(e[i, ])
      check(is.na(expected_parent[child]), "atom has multiple earlier tree neighbours")
      expected_parent[child] <- par
    }
    check(identical(x$parent, expected_parent), "parent does not match rooted ordered tree")
    check(identical(x$children, lapply(seq_len(k), function(i) which(expected_parent == i))), "children do not match parents")
    nonroots <- which(!is.na(expected_parent))
    check(identical(x$separator_atoms, nonroots), "separator-to-atom mapping mismatch")
    check(length(x$separators_index) == length(nonroots), "ordered separator count mismatch")
    expected_ss <- lapply(nonroots, function(i) intersect(aa[[i]], aa[[expected_parent[i]]]))
    check(identical(x$separators_index, expected_ss), "separator differs from child-parent intersection")
    check(all(vapply(expected_ss, function(s) length(s) > 0L && .is_complete(g, s), logical(1))), "separator is empty or not a clique")
    expected_rr <- lapply(seq_len(k), function(i) if (is.na(expected_parent[i])) aa[[i]] else setdiff(aa[[i]], aa[[expected_parent[i]]]))
    check(identical(x$residuals_index, expected_rr), "residual definitions mismatch")
    check(identical(x$complete, vapply(aa, function(v) .is_complete(g, v), logical(1))), "complete flags mismatch")
    oo <- x$component_orderings_index
    check(length(oo) == nc, "component ordering count mismatch")
    check(identical(x$component_orderings, lapply(oo, .external_sets, norm = norm)), "external ordering labels mismatch")
    check(identical(x$prime_ordering, if (nc == 1L) x$component_orderings[[1L]] else if (!nc) list() else NULL), "prime_ordering component representation mismatch")
    for (ci in seq_len(nc)) {
      ids <- which(ac == ci); ordering <- oo[[ci]]
      check(length(ordering) == 2L * length(ids) - 1L, paste("component", ci, "perfect ordering length mismatch"))
      previous <- integer()
      for (j in seq_along(ids)) {
        i <- ids[j]; a <- aa[[i]]
        check(identical(ordering[[2L * j - 1L]], a), paste("atom", i, "alternating ordering mismatch"))
        if (j > 1L) {
          s <- ordering[[2L * j - 2L]]
          check(setequal(s, intersect(a, previous)), paste("atom", i, "separator differs from previous-union intersection"))
          par <- expected_parent[i]
          check(!is.na(par) && par < i && par %in% ids[seq_len(j - 1L)], paste("atom", i, "has no valid earlier parent"))
          if (!is.na(par)) check(setequal(s, intersect(a, aa[[par]])) && all(s %in% aa[[par]]), paste("atom", i, "ordering separator differs from parent intersection"))
        } else check(is.na(expected_parent[i]), "component first atom is not a root")
        previous <- union(previous, a)
      }
    }
    # In a forest, the induced holder graph is connected iff it has h-1 edges.
    holders <- vector("list", p)
    for (i in seq_along(aa)) for (v in aa[[i]]) holders[[v]] <- c(holders[[v]], i)
    for (v in seq_len(p)) {
      ids <- holders[[v]]
      check(sum(e[, 1L] %in% ids & e[, 2L] %in% ids) == length(ids) - 1L, paste("running intersection fails at vertex", v))
    }
    if (expensive) {
      ge <- igraph::as_edgelist(g, names = FALSE)
      for (i in seq_len(nrow(ge))) check(any(vapply(aa, function(a) all(ge[i, ] %in% a), logical(1))), paste("uncovered graph edge", i))
      for (i in seq_len(k)) for (j in seq_len(k)) if (i != j)
        check(!all(aa[[i]] %in% aa[[j]]), paste("atom", i, "is contained in atom", j))
    }
  }, error = function(e) errors <<- c(errors, paste("Malformed decomposition:", conditionMessage(e))))
  if (!length(errors)) TRUE else structure(list(valid = FALSE, errors = unique(errors)), class = "atom_decomp_validation")
}
