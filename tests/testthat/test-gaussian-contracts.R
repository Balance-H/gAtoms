test_that("Both public Gaussian functions reject every nonfinite value before fitting", {
  local_mocked_bindings(.fit_gaussian_atom = function(...) stop("unexpected fitting"))
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    X <- matrix(seq_len(20), 5, 4); X[1,1] <- bad
    expect_error(gaussian_suffstats(X), "only finite, non-missing")
    expect_error(fit_gaussian_graph(X, igraph::make_ring(4)), "only finite, non-missing")
  }
  expect_named(formals(gaussian_suffstats), c("data", "mean"))
  expect_named(formals(fit_gaussian_graph), c("data","graph","decomp","engine","covariance","mean","control"))
})

for (backend in c("gRim", "ggm")) local({
  engine <- backend
  test_that(paste(engine, "preserves its covariance with the required inversion count"), {
    skip_if_not_installed(engine)
    Sigma <- diag(4); Sigma[1,3] <- Sigma[3,1] <- 1e-12
    K <- chol2inv(chol(Sigma))
    if (engine == "gRim") {
      local_mocked_bindings(fit_ggm_grips = function(...) list(K=K, Sigma=Sigma,
        converged=TRUE, details=list(iter=2L)), .package="gRim")
    } else {
      local_mocked_bindings(fitConGraph = function(...) list(Shat=Sigma,it=2L), .package="ggm")
    }
    inversions <- 0L
    env <- new.env(parent=environment(.fit_gaussian_atom))
    env$chol2inv <- function(...) { inversions <<- inversions+1L; base::chol2inv(...) }
    atom <- .fit_gaussian_atom; environment(atom) <- env
    adapter <- if (engine == "gRim") .fit_gaussian_grim else .fit_gaussian_ggm
    environment(adapter) <- env
    assign(if (engine == "gRim") ".fit_gaussian_grim" else ".fit_gaussian_ggm", adapter, env)
    graph <- igraph::make_ring(4); control <- .gaussian_control(list())
    raw <- adapter(diag(4), graph, 100, control)
    expect_identical(raw$Sigma, Sigma)
    expect_equal(raw$K, K)
    expect_identical(inversions, if (engine == "gRim") 0L else 1L)
    inversions <- 0L
    fit <- atom(diag(4), graph, 100, FALSE, engine, control)
    expect_identical(inversions, if (engine == "gRim") 0L else 1L)
    expect_identical(fit$Sigma, Sigma)
    expect_identical(fit$K[1,3], 0)
    expect_equal(fit$K %*% fit$Sigma, diag(4), tolerance=control$contract_tol)
    adj <- as.matrix(igraph::as_adjacency_matrix(graph)) != 0; diag(adj) <- TRUE
    expect_equal(fit$Sigma[adj], diag(4)[adj])
    inversions <- 0L
    atom(diag(4), igraph::make_full_graph(4), 100, TRUE, engine, control)
    expect_identical(inversions, 1L)
    inversions <- 0L
    separator <- .fit_gaussian_separator; environment(separator) <- env
    separator(diag(2), 1:2, 1L)
    expect_identical(inversions, 1L)
  })
})

test_that("Local matrix contracts reject invalid covariance and material graph residuals", {
  raw <- list(K=diag(4),Sigma=diag(4),converged=TRUE,iterations=1L,method="mock",message="")
  local_mocked_bindings(.fit_gaussian_grim = function(...) raw)
  run <- function() .fit_gaussian_atom(diag(4),igraph::make_ring(4),100,FALSE,"gRim",.gaussian_control(list()))
  for (bad in list(matrix(1,2,2),diag(c(NA,1,1,1)),diag(c(Inf,1,1,1)),diag(c(-1,1,1,1)))) {
    raw$Sigma <- bad
    expect_error(run(),class="atom_gaussian_error")
  }
  raw$Sigma <- diag(4); raw$Sigma[1,2] <- 0.1
  expect_error(run(),"symmetry")
  raw$Sigma <- diag(4)*2
  expect_error(run(),"numerically inconsistent")
  raw$K <- diag(4)/2
  expect_error(run(),"score equations")
  raw$Sigma <- diag(c(1e-8,1,1,1)); raw$K <- diag(c(1e8,1,1,1))
  raw$K[2,4] <- raw$K[4,2] <- 1e-4
  e <- tryCatch(run(),atom_gaussian_error=identity)
  expect_s3_class(e,"atom_gaussian_error")
  expect_identical(e$status,"nonconverged")
  expect_match(conditionMessage(e),"graph zeros")
})

test_that("An explicit engine nonconvergence flag retains its status", {
  local_mocked_bindings(.fit_gaussian_grim = function(...) list(converged=FALSE))
  e <- tryCatch(.fit_gaussian_atom(diag(4),igraph::make_ring(4),100,FALSE,
    "gRim",.gaussian_control(list())),atom_gaussian_error=identity)
  expect_s3_class(e,"atom_gaussian_error")
  expect_identical(e$status,"nonconverged")
})
