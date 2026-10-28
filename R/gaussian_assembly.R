.assemble_gaussian_precision <- function(local_fits, separator_fits, p, labels) {
  blocks <- c(local_fits, separator_fits)
  signs <- c(rep(1, length(local_fits)), rep(-1, length(separator_fits)))
  entries <- lapply(seq_along(blocks), function(i) {
    z <- blocks[[i]]; nz <- which(upper.tri(z$K, diag = TRUE) & z$K != 0, arr.ind = TRUE)
    list(i = z$vertices[nz[, 1L]], j = z$vertices[nz[, 2L]], x = signs[i]*z$K[nz])
  })
  Matrix::drop0(Matrix::sparseMatrix(i = unlist(lapply(entries, `[[`, "i")),
    j = unlist(lapply(entries, `[[`, "j")), x = unlist(lapply(entries, `[[`, "x")),
    dims = c(p, p), dimnames = list(labels, labels), symmetric = TRUE))
}

.check_gaussian_separators <- function(decomp, local_fits, separator_fits, tol) {
  for (j in seq_along(separator_fits)) {
    z <- separator_fits[[j]]; child <- decomp$separator_atoms[j]
    for (i in c(child, decomp$parent[child])) {
      pos <- match(z$vertices, local_fits[[i]]$vertices)
      block <- local_fits[[i]]$Sigma[pos, pos, drop = FALSE]
      scale <- sqrt(outer(diag(z$Sigma), diag(z$Sigma)))
      if (max(abs(block-z$Sigma)/scale) > tol)
        .gaussian_failure("solver_failed", "Local separator covariance violates the canonical covariance contract.", i, z$vertices)
    }
  }
  invisible(TRUE)
}

.assemble_gaussian_covariance <- function(decomp, local_fits, separator_fits,
                                          precision, mode) {
  if (mode == "none") return(NULL)
  if (mode == "solve") return(as.matrix(Matrix::solve(precision)))
  p <- nrow(precision)
  Sigma <- matrix(0, p, p, dimnames = dimnames(precision))
  for (ci in unique(decomp$components$atoms)) {
    history <- integer()
    for (i in which(decomp$components$atoms == ci)) {
      fit <- local_fits[[i]]; v <- fit$vertices
      if (is.na(decomp$parent[i])) {
        Sigma[v, v] <- fit$Sigma
      } else {
        sep <- separator_fits[[match(i, decomp$separator_atoms)]]
        s <- sep$vertices; r <- decomp$residuals_index[[i]]; b <- setdiff(history, s)
        ri <- match(r, v); si <- match(s, v)
        Sigma[r, r] <- fit$Sigma[ri, ri, drop = FALSE]
        Sigma[r, s] <- fit$Sigma[ri, si, drop = FALSE]
        Sigma[s, r] <- t(Sigma[r, s, drop = FALSE])
        Sigma[s, s] <- sep$Sigma
        if (length(b)) {
          cross <- fit$Sigma[ri, si, drop = FALSE] %*% sep$K %*% Sigma[s, b, drop = FALSE]
          Sigma[r, b] <- cross; Sigma[b, r] <- t(cross)
        }
      }
      history <- union(history, v)
    }
  }
  Sigma
}
