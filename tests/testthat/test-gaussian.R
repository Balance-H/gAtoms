gaussian_example <- function(p, n = 250L) {
  set.seed(621)
  X <- matrix(rnorm(n*p), n, p)
  X <- X + matrix(rnorm(n), n, p)
  sweep(X, 2L, seq_len(p)/3, "+")
}
gaussian_graph <- function(kind) {
  switch(kind,
    path = igraph::make_ring(5, circular = FALSE),
    prime = igraph::make_ring(4),
    mixed = igraph::add_edges(igraph::add_vertices(igraph::make_ring(4), 2), c(3,5,4,5,5,6)))
}
gaussian_direct <- function(X, G, engine) {
  ss <- gaussian_suffstats(X)
  z <- .fit_gaussian_atom(ss$S_ml, G, ss$n, FALSE, engine, .gaussian_control(list()))
  list(K = z$K, Sigma = z$Sigma,
    ll = .gaussian_loglik(ss$n, ss$p, -2*sum(log(diag(chol(z$Sigma)))), sum(ss$S_ml*t(z$K))))
}

test_that("Gaussian statistics use one global SSD and the ML divisor", {
  X <- gaussian_example(3)
  z <- gaussian_suffstats(X)
  expect_equal(z$ssd, crossprod(sweep(X, 2, colMeans(X))))
  expect_equal(z$S_ml, stats::cov(X)*(nrow(X)-1)/nrow(X))
  expect_equal(z$S_unbiased, stats::cov(X))
  z0 <- gaussian_suffstats(X, mean = "zero")
  expect_equal(z0$S_ml, crossprod(X)/nrow(X))
  expect_equal(z0$mean, rep(0, 3))
  expect_equal(z0$S_unbiased, z0$S_ml)
  expect_named(z, c("n","p","mean","ssd","S_ml","S_unbiased","mean_mode"))
  expect_error(gaussian_suffstats(data.frame(x = letters)), "numeric")
})

test_that("Complete atoms use closed form and likelihood constants are correct", {
  X <- gaussian_example(5); G <- igraph::make_full_graph(5)
  for (mu in c("estimate", "zero")) {
    f <- fit_gaussian_graph(X, G, mean = mu, engine = "ggm")
    S <- f$suffstats$S_ml
    expect_s3_class(f, "atom_gaussian_fit")
    expect_s4_class(precision(f), "symmetricMatrix")
    expect_equal(unname(vcov(f)), unname(S), tolerance = 1e-12)
    expect_equal(unname(as.matrix(precision(f))), unname(solve(S)), tolerance = 1e-10)
    expect_identical(f$local_fits[[1]]$engine, "closed_form")
    expect_identical(f$local_fits[[1]]$iterations, 0L)
    ll <- -nrow(X)/2 * (5*log(2*pi) + as.numeric(determinant(S, logarithm=TRUE)$modulus) + 5)
    expect_equal(as.numeric(logLik(f)), ll, tolerance = 1e-10)
    df <- 15 + if (mu == "estimate") 5 else 0
    expect_equal(attr(logLik(f), "df"), df)
    expect_equal(AIC(f), -2*ll+2*df)
    expect_equal(BIC(f), -2*ll+log(nrow(X))*df)
    expect_equal(unname(f$h), as.numeric(precision(f) %*% f$mean))
    expect_equal(as.numeric(deviance(f)), 0)
    expect_output(print(f), "Gaussian graph fit")
    expect_equal(summary(f)$n, nrow(X))
  }
})

test_that("Ordered separator assembly, recursive covariance and scalar summaries agree", {
  G <- gaussian_graph("path"); X <- gaussian_example(5)
  d <- atom_decompose(G)
  f <- fit_gaussian_graph(X, G, decomp = d)
  expect_identical(f$decomp, d)
  expect_equal(vcov(f), as.matrix(Matrix::solve(precision(f))), tolerance = 1e-10)
  expect_equal(f$logdet$precision, as.numeric(determinant(as.matrix(precision(f)), logarithm=TRUE)$modulus), tolerance=1e-10)
  expect_equal(f$likelihood$trace, sum(f$suffstats$S_ml*t(as.matrix(precision(f)))), tolerance=1e-10)
  expect_equal(f$likelihood$trace, 5, tolerance=1e-8)
  expect_equal(attr(deviance(f), "df"), choose(5,2)-4)
  expect_equal(fit_gaussian_graph(X,G,covariance="solve")$covariance, vcov(f), tolerance=1e-10)
  none <- fit_gaussian_graph(X,G,covariance="none")
  expect_null(vcov(none))
  expect_equal(as.matrix(precision(none)), as.matrix(precision(f)))
  expect_equal(logLik(none), logLik(f))
  bad <- f$local_fits; bad[[d$separator_atoms[1]]]$Sigma <- bad[[d$separator_atoms[1]]]$Sigma * 2
  expect_error(.check_gaussian_separators(d,bad,f$separator_fits,1e-6), class="atom_gaussian_error")
})

test_that("Disconnected components and names are preserved", {
  G <- igraph::disjoint_union(igraph::make_ring(3,circular=FALSE), igraph::make_empty_graph(1,directed=FALSE))
  G <- igraph::set_vertex_attr(G,"name",value=c("a b","x","z","last"))
  X <- gaussian_example(4); colnames(X) <- igraph::V(G)$name
  f <- fit_gaussian_graph(X[,4:1],G)
  ref <- fit_gaussian_graph(X,G)
  expect_equal(precision(f),precision(ref))
  expect_equal(f$suffstats$n,nrow(X))
  expect_equal(colnames(precision(f)),igraph::V(G)$name)
  expect_equal(as.numeric(precision(f)[1:3,4]),rep(0,3))
  expect_equal(as.numeric(vcov(f)[1:3,4]),rep(0,3))
  expect_equal(vcov(f),as.matrix(Matrix::solve(precision(f))),tolerance=1e-10)
  one <- fit_gaussian_graph(matrix(1:10,ncol=1),igraph::make_empty_graph(1,directed=FALSE))
  expect_equal(dim(vcov(one)),c(1L,1L))
  expect_error(fit_gaussian_graph(X,G,decomp=atom_decompose(igraph::make_full_graph(4))), "same graph")
  bad <- atom_decompose(G); bad$parent[1] <- 1L
  expect_error(fit_gaussian_graph(X,G,decomp=bad), "Invalid decomp")
  expect_error(fit_gaussian_graph(X,G,control=list(ridge=1e-3)),"Unknown")
  expect_error(fit_gaussian_graph(X,G,control=list(tol=NA)),"positive")
})

test_that("Singular statistics fail structurally without silent ridge", {
  X <- cbind(1:10,1:10)
  e <- tryCatch(fit_gaussian_graph(X,igraph::make_full_graph(2)),atom_gaussian_error=identity)
  expect_s3_class(e,"atom_gaussian_error")
  expect_identical(e$status,"singular_statistics")
  expect_equal(e$atom_id,1L)
  expect_error(fit_gaussian_graph(matrix(1,1,3),igraph::make_full_graph(3)),class="atom_gaussian_error")
  expect_error(fit_gaussian_graph(gaussian_example(3,n=3),igraph::make_full_graph(3)),class="atom_gaussian_error")
  expect_error(fit_gaussian_graph(gaussian_example(3,n=2),igraph::make_full_graph(3),mean="zero"),class="atom_gaussian_error")
  # A constrained MLE can exist while the saturated reference does not.
  X <- matrix(c(1,2,3,2,4,1,3,1,4),3,3)
  X <- cbind(X,X[,1])
  f <- fit_gaussian_graph(X,igraph::make_empty_graph(4,directed=FALSE))
  expect_true(is.na(deviance(f)))
  expect_match(f$diagnostics$deviance,"singular")
})

for (eng in c("gRim","ggm")) local({
  engine <- eng
  test_that(paste(engine,"matches full-graph fits for decomposable, mixed and prime graphs"), {
    skip_if_not_installed(engine)
    for (kind in c("path","mixed","prime")) {
      G <- gaussian_graph(kind); X <- gaussian_example(igraph::vcount(G))
      f <- fit_gaussian_graph(X,G,engine=engine)
      ref <- gaussian_direct(X,G,engine)
      if (kind == "mixed") {
        expect_gt(length(f$local_fits),1L)
        expect_true(any(!f$decomp$complete))
        expect_true(any(vapply(f$local_fits,function(z) z$engine==engine,logical(1))))
      }
      if (kind == "prime") expect_length(f$local_fits,1L)
      expect_equal(unname(as.matrix(precision(f))),unname(ref$K),tolerance=1e-6)
      expect_equal(unname(vcov(f)),unname(ref$Sigma),tolerance=1e-6)
      expect_equal(as.numeric(logLik(f)),ref$ll,tolerance=1e-7)
      expect_equal(vcov(f),as.matrix(Matrix::solve(precision(f))),tolerance=1e-6)
      expect_equal(f$logdet$precision,as.numeric(determinant(as.matrix(precision(f)))$modulus),tolerance=1e-6)
      expect_equal(f$likelihood$trace,sum(f$suffstats$S_ml*t(as.matrix(precision(f)))),tolerance=1e-8)
      expect_equal(f$likelihood$trace,igraph::vcount(G),tolerance=1e-6)
    }
  })
})

test_that("Both engines retain the supplied ML covariance scale", {
  skip_if_not_installed("gRim"); skip_if_not_installed("ggm")
  G <- gaussian_graph("mixed"); X <- gaussian_example(6)
  a <- fit_gaussian_graph(X,G,engine="gRim")
  b <- fit_gaussian_graph(X,G,engine="ggm")
  expect_equal(precision(a),precision(b),tolerance=1e-6)
  expect_equal(vcov(a),vcov(b),tolerance=1e-6)
  S <- gaussian_suffstats(X)$S_ml
  dimnames(S) <- list(paste0("v",1:6),paste0("v",1:6))
  adj <- as.matrix(igraph::as_adjacency_matrix(G)); dimnames(adj) <- dimnames(S)
  raw <- ggm::fitConGraph(adj,S,nrow(X),tol=1e-10)
  expect_equal(unname(raw$Shat),unname(vcov(b)),tolerance=1e-6)
  raw2 <- gRim::fit_ggm_grips(S,t(igraph::as_edgelist(G,names=FALSE)),nobs=nrow(X),eps=1e-10)
  expect_equal(unname(raw2$Sigma),unname(vcov(a)),tolerance=1e-6)
})

test_that("Iteration-limited inconsistent engine outputs are rejected", {
  skip_if_not_installed("gRim")
  G <- gaussian_graph("prime"); X <- gaussian_example(4)
  e <- tryCatch(fit_gaussian_graph(X,G,control=list(maxit=1L,tol=1e-15,fallback=FALSE)),atom_gaussian_error=identity)
  expect_s3_class(e,"atom_gaussian_error")
  expect_identical(e$status,"solver_failed")
  expect_match(conditionMessage(e),"numerically inconsistent")
})

test_that("Failed local engines never produce partial global matrices", {
  local_mocked_bindings(.fit_gaussian_grim = function(...) stop("simulated solver failure"),
    .assemble_gaussian_precision = function(...) stop("assembly must not run"))
  e <- tryCatch(fit_gaussian_graph(gaussian_example(6),gaussian_graph("mixed")),atom_gaussian_error=identity)
  expect_s3_class(e,"atom_gaussian_error")
  expect_identical(e$status,"solver_failed")
  expect_match(conditionMessage(e),"simulated solver failure")
  expect_true(length(e$vertices) > 0)
})

test_that("The default engine is gRim and timing$decomposition tracks atom_decompose only", {
  G <- gaussian_graph("path"); X <- gaussian_example(igraph::vcount(G))
  dec <- atom_decompose(G)

  # default is equivalent to engine = "gRim"
  a <- fit_gaussian_graph(X, G, decomp = dec)
  b <- fit_gaussian_graph(X, G, decomp = dec, engine = "gRim")
  expect_identical(precision(a), precision(b))
  expect_identical(vcov(a), vcov(b))
  expect_identical(as.numeric(logLik(a)), as.numeric(logLik(b)))
  expect_true(all(vapply(a$local_fits, function(z) z$engine != "auto", logical(1))))

  # reusing a supplied decomposition costs zero decomposition time
  expect_identical(a$timing$decomposition, 0)

  # with decomp = NULL only atom_decompose() is timed
  c0 <- fit_gaussian_graph(X, G)
  expect_true(is.numeric(c0$timing$decomposition))
  expect_length(c0$timing$decomposition, 1L)
  expect_true(is.finite(c0$timing$decomposition))
  expect_true(c0$timing$decomposition >= 0)
  expect_identical(as.numeric(logLik(c0)), as.numeric(logLik(a)))
  expect_identical(c0$decomp$atoms_index, dec$atoms_index)

  # engine matching is validated
  expect_error(fit_gaussian_graph(X, G, engine = "auto"), "arg")
})

test_that("Complete graphs never dispatch to optional engines", {
  local_mocked_bindings(.fit_gaussian_grim = function(...) stop("unexpected engine"),
    .fit_gaussian_ggm = function(...) stop("unexpected engine"))
  for (engine in c("gRim","ggm")) {
    f <- fit_gaussian_graph(gaussian_example(6),igraph::make_full_graph(6),engine=engine)
    expect_identical(f$local_fits[[1]]$status,"closed_form")
  }
})
