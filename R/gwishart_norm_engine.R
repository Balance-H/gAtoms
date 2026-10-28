# gAtoms G-Wishart normalizing constant -- local Monte Carlo engine adapter.
#
# This file contains the ONLY place where BDgraph is used.  It is responsible
# for the local, incomplete prime graphs only; everything else in the module is
# exact.  BDgraph::gnorm() returns an approximation of log I_{G_P}(b, D_P) for
# the parameterization documented in gwishart_norm_statistics.R, so any local
# result produced here is flagged approximate = TRUE.

# BDgraph is an optional runtime dependency (DESCRIPTION is intentionally not
# modified).  Only a graph with at least one incomplete prime needs it.
.gwishart_bdgraph_available <- function() {
  requireNamespace("BDgraph", quietly = TRUE) &&
    exists("gnorm", asNamespace("BDgraph"), inherits = FALSE)
}

# Build the strict upper-triangular binary adjacency matrix expected by
# BDgraph::gnorm(): zero diagonal and a zero lower triangle.  The induced
# subgraph is used so that the local vertex order is exactly the order of the
# corresponding principal block of D.
.gwishart_bdgraph_adjacency <- function(graph, index) {
  adj <- .gwishart_adjacency(graph, index)
  q <- nrow(adj)
  if (!identical(dim(adj), c(q, q)) || !all(adj %in% c(0, 1)))
    stop("Local adjacency must be a square 0/1 matrix.", call. = FALSE)
  adj[lower.tri(adj, diag = TRUE)] <- 0
  if (any(diag(adj) != 0) || any(adj[lower.tri(adj)] != 0))
    stop("Local adjacency must be strict upper triangular.", call. = FALSE)
  adj
}

# Call BDgraph::gnorm() on one induced prime graph and validate its output.
# Every failure is reported with the local context required for diagnosis.
.gwishart_bdgraph_log_norm <- function(graph, index, D, b, iter,
                                       atom_id = NA_integer_,
                                       vertices = names(D)) {
  start <- proc.time()[[3L]]
  if (!.gwishart_bdgraph_available())
    stop(sprintf(paste0("BDgraph is required for incomplete prime atom %s ",
      "(local vertices: %s; local dimension: %d; b = %s; iter = %d) but the ",
      "optional package BDgraph (or BDgraph::gnorm) is not available."),
      format(atom_id), paste(vertices, collapse = ","), ncol(D), format(b),
      as.integer(iter)), call. = FALSE)
  adj <- .gwishart_bdgraph_adjacency(graph, index)
  if (!identical(dim(adj), dim(D)))
    stop(sprintf(paste0("Internal ordering mismatch for atom %s (local ",
      "vertices: %s): adjacency dimension %d differs from D dimension %d."),
      format(atom_id), paste(vertices, collapse = ","), nrow(adj), nrow(D)),
      call. = FALSE)
  warnings <- character()
  value <- tryCatch(
    withCallingHandlers(
      BDgraph::gnorm(adj, b = as.numeric(b), D = D, iter = as.integer(iter)),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
    error = function(e) {
      stop(sprintf(paste0("BDgraph::gnorm failed for incomplete prime atom %s ",
        "(local vertices: %s; local dimension: %d; b = %s; iter = %d): %s"),
        format(atom_id), paste(vertices, collapse = ","), ncol(D), format(b),
        as.integer(iter), conditionMessage(e)), call. = FALSE)
    })
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || is.nan(value) ||
      !is.finite(value))
    stop(sprintf(paste0("BDgraph::gnorm returned an invalid log normalizing ",
      "constant for incomplete prime atom %s (local vertices: %s; local ",
      "dimension: %d; b = %s; iter = %d): got %s."),
      format(atom_id), paste(vertices, collapse = ","), ncol(D), format(b),
      as.integer(iter), paste(format(value), collapse = ", ")), call. = FALSE)
  list(log_norm = as.numeric(value), adjacency = adj,
       iterations = as.integer(iter), message = paste(warnings, collapse = "; "),
       elapsed = unname(proc.time()[[3L]] - start), status = "converged",
       engine = "BDgraph", approximate = TRUE)
}

# Normalized local result object for one prime component.  The stored
# `adjacency` is always the strict upper-triangular binary form used for the
# BDgraph contract, so every local fit exposes the same representation.
.gwishart_local_fit <- function(atom_id, component_id, index, graph, vertices,
                                D, b, complete, iter) {
  start <- proc.time()[[3L]]
  local_D <- .gwishart_block(D, index)
  if (complete) {
    list(atom_id = as.integer(atom_id), component_id = as.integer(component_id),
         vertices = vertices, size = length(index), complete = TRUE,
         engine = "closed_form", approximate = FALSE,
         log_norm = .gwishart_log_norm_complete(local_D, b),
         elapsed = unname(proc.time()[[3L]] - start),
         status = "closed_form", message = "",
         adjacency = .gwishart_bdgraph_adjacency(graph, index))
  } else {
    z <- .gwishart_bdgraph_log_norm(graph, index, local_D, b, iter,
                                   atom_id = atom_id, vertices = vertices)
    list(atom_id = as.integer(atom_id), component_id = as.integer(component_id),
         vertices = vertices, size = length(index), complete = FALSE,
         engine = "BDgraph", approximate = TRUE, log_norm = z$log_norm,
         elapsed = z$elapsed, status = z$status, message = z$message,
         adjacency = z$adjacency)
  }
}

# Normalized local result object for one clique separator.  Separators of a
# prime decomposition are complete by construction, so this is always exact.
.gwishart_separator_fit <- function(separator_id, component_id, index, vertices, D, b) {
  start <- proc.time()[[3L]]
  list(separator_id = as.integer(separator_id),
       component_id = as.integer(component_id), vertices = vertices,
       size = length(index), engine = "closed_form", approximate = FALSE,
       log_norm = .gwishart_log_norm_complete(.gwishart_block(D, index), b),
       elapsed = unname(proc.time()[[3L]] - start), status = "closed_form",
       message = "")
}
