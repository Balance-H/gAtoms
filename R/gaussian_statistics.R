#' Canonical Gaussian sufficient statistics
#' @param data A numeric matrix or all-numeric data frame; rows are observations.
#' @param mean Estimate the mean, or fix it at zero.
#' @return A list with n, p, mean, ssd, S_ml, S_unbiased, and mean_mode.
#' @details S_ml is always SSD/n. With estimated mean, S_unbiased is SSD/(n-1)
#'   (NA when n=1). With known zero mean both covariance estimates are SSD/n.
#'   The Gaussian module currently requires complete finite observations.
#' @export
gaussian_suffstats <- function(data, mean = c("estimate", "zero")) {
  mean <- match.arg(mean)
  if (is.data.frame(data)) {
    if (!all(vapply(data, is.numeric, logical(1)))) stop("All data columns must be numeric.", call. = FALSE)
    data <- as.matrix(data)
  }
  if (!is.matrix(data) || !is.numeric(data) || !ncol(data))
    stop("data must be a numeric matrix with at least one column.", call. = FALSE)
  if (any(!is.finite(data))) stop("data must contain only finite, non-missing values.", call. = FALSE)
  n <- nrow(data); p <- ncol(data)
  if (!n) stop("No observations remain.", call. = FALSE)
  mu <- if (mean == "estimate") colMeans(data) else rep(0, p)
  names(mu) <- colnames(data)
  centered <- if (mean == "estimate") sweep(data, 2L, mu) else data
  ssd <- crossprod(centered)
  if (any(!is.finite(ssd))) stop("Sufficient statistics overflowed; rescale data.", call. = FALSE)
  list(n = n, p = p, mean = mu, ssd = ssd, S_ml = ssd/n,
       S_unbiased = if (mean == "zero") ssd/n else if (n > 1L) ssd/(n-1L) else ssd * NA_real_,
       mean_mode = mean)
}

.gaussian_failure <- function(status, message, atom_id = NA_integer_, vertices = integer()) {
  stop(structure(list(message = message, call = NULL, status = status,
    atom_id = atom_id, vertices = vertices),
    class = c("atom_gaussian_error", "error", "condition")))
}

.gaussian_chol <- function(S, status = "singular_statistics") {
  z <- tryCatch(chol(S), error = function(e) NULL)
  if (is.null(z) || any(!is.finite(z)) || any(diag(z) <= 0))
    .gaussian_failure(status, "Required covariance is not positive definite; no regularization was applied.")
  z
}

.gaussian_logdet <- function(local_fits, separator_fits) {
  ld <- function(z) 2 * sum(log(diag(.gaussian_chol(z$Sigma, "solver_failed"))))
  s <- sum(vapply(local_fits, ld, 0.0)) - sum(vapply(separator_fits, ld, 0.0))
  list(precision = -s, covariance = s)
}

.gaussian_trace <- function(suffstats, local_fits, separator_fits) {
  tr <- function(z) sum(suffstats$S_ml[z$vertices, z$vertices, drop = FALSE] * t(z$K))
  sum(vapply(local_fits, tr, 0.0)) - sum(vapply(separator_fits, tr, 0.0))
}

.gaussian_loglik <- function(n, p, logdet, trace) n/2 * (logdet - trace - p * log(2*pi))
