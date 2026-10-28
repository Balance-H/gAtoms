ll_graph <- function(kind) {
  g <- switch(kind,path=igraph::make_ring(4,circular=FALSE),prime=igraph::make_ring(4),
    mixed=igraph::add_edges(igraph::add_vertices(igraph::make_ring(4),2),c(3,5,4,5,5,6)),
    complete=igraph::make_full_graph(3),disconnected=igraph::disjoint_union(
      igraph::make_ring(3,circular=FALSE),igraph::make_empty_graph(1,directed=FALSE)))
  igraph::set_vertex_attr(g,'name',value=LETTERS[seq_len(igraph::vcount(g))])
}
ll_data <- function(g) {
  grid <- do.call(expand.grid,stats::setNames(rep(list(c('lo','hi')),igraph::vcount(g)),igraph::V(g)$name))
  grid[rep(seq_len(nrow(grid)),1L+(seq_len(nrow(grid))*7L)%%9L),,drop=FALSE]
}
ll_oracle <- function(dat,g) {
  tab <- do.call(table,dat)
  generators <- lapply(igraph::max_cliques(g),as.integer)
  z <- stats::loglin(tab,generators,fit=TRUE,eps=1e-10,iter=10000,print=FALSE)
  positive <- tab>0
  ll <- sum(tab[positive]*log(z$fit[positive]/sum(tab)))
  sat <- sum(tab[positive]*log(tab[positive]/sum(tab)))
  list(table=tab,prob=z$fit/sum(tab),ll=ll,dev=2*(sat-ll),generators=generators)
}

test_that('Complete atoms and separators bypass iterative engines', {
  local_mocked_bindings(.fit_loglinear_loglin=function(...) stop('unexpected engine'),
    .fit_loglinear_grim=function(...) stop('unexpected engine'))
  for (kind in c('complete','path')) {
    g <- ll_graph(kind); dat <- ll_data(g); f <- fit_loglinear_graph(dat,g)
    expect_true(all(f$decomp$complete))
    expect_true(all(vapply(f$local_fits,function(z) z$engine=='closed_form' && z$iterations==0,logical(1))))
    for (z in c(f$local_fits,f$separator_fits)) expect_equal(z$fitted_counts,z$observed)
    if (kind=='complete') {
      expect_length(f$local_fits,1)
      expect_equal(as.numeric(logLik(f)),f$likelihood$saturated)
      expect_equal(deviance(f),0)
    }
  }
})

for (eng in c('loglin','gRim')) local({
  engine <- eng
  test_that(paste(engine,'matches full-table graphical fits'), {
    if (engine=='gRim') skip_if_not_installed('gRim')
    for (kind in c('path','mixed','prime','disconnected')) {
      g <- ll_graph(kind); dat <- ll_data(g); dec <- atom_decompose(g)
      f <- fit_loglinear_graph(dat[,rev(names(dat))],g,decomp=dec,engine=engine)
      ref <- ll_oracle(dat,g)
      expect_identical(f$decomp,dec)
      expect_equal(as.numeric(logLik(f)),ref$ll,tolerance=1e-6)
      expect_equal(deviance(f),ref$dev,tolerance=1e-6)
      full <- fitted(f,type='full')
      expect_equal(as.vector(full),as.vector(ref$prob),tolerance=1e-6)
      expect_equal(sum(full),1,tolerance=1e-6)
      for (s in ref$generators) expect_equal(as.vector(margin.table(full,s)),
        as.vector(margin.table(ref$table,s)/sum(ref$table)),tolerance=1e-6)
      grid <- do.call(expand.grid,f$levels)
      expect_equal(as.vector(full),exp(.loglinear_logprob(fitted(f),grid)),tolerance=1e-12)
      expect_equal(sum(.loglinear_logprob(fitted(f),dat)),as.numeric(logLik(f)),tolerance=1e-6)
      if (kind=='mixed') {
        expect_gt(length(f$local_fits),1)
        expect_true(any(!dec$complete))
        expect_true(any(vapply(f$local_fits,function(z) z$engine==engine,logical(1))))
      }
      if (kind=='prime') expect_length(f$local_fits,1)
      expect_output(print(f),'Graphical log-linear fit')
      expect_equal(summary(f)$n,nrow(dat))
    }
  })
})

test_that('Raw and sufficient-margin pathways agree and validate compatibility', {
  g <- ll_graph('mixed'); dat <- ll_data(g); dec <- atom_decompose(g)
  ss <- margin_counts(dat,c(atoms(dec),separators(dec)))
  a <- fit_loglinear_graph(dat,g,engine='loglin')
  b <- fit_loglinear_graph(ss,g,decomp=dec,engine='loglin')
  expect_equal(logLik(a),logLik(b))
  expect_equal(deviance(a),deviance(b))
  expect_equal(fitted(a),fitted(b))
  ss$observed_counts <- NULL
  b <- fit_loglinear_graph(ss,g,engine='loglin')
  expect_true(is.na(deviance(b)))
  expect_match(b$diagnostics$deviance,'unavailable')
  expect_error(fit_loglinear_graph(margin_counts(dat,list('A')),g),'lack')
  bad <- ss; bad$margins[[1]][1] <- -1
  expect_error(fit_loglinear_graph(bad,g),class='atom_loglinear_error')
  expect_error(fit_loglinear_graph(dat,g,decomp=atom_decompose(ll_graph('complete'))),'match')
  names(dat)[1] <- 'unknown'
  expect_error(fit_loglinear_graph(dat,g),'match')
})

test_that('Sampling zeros, unused levels and separator zero support are preserved', {
  for (engine in c('loglin',if (requireNamespace('gRim',quietly=TRUE)) 'gRim')) for (kind in c('path','prime')) {
    g <- ll_graph(kind); dat <- ll_data(g)
    dat <- dat[dat$A=='lo' | dat$B=='hi',,drop=FALSE]
    dat$A <- factor(dat$A,levels=c('hi','unused','lo'))
    f <- fit_loglinear_graph(dat,g,engine=engine,control=list(maxit=100L))
    expect_identical(f$levels$A,c('hi','unused','lo'))
    expect_equal(sum(fitted(f,type='full')),1,tolerance=1e-6)
    expect_true(all(is.finite(.loglinear_logprob(fitted(f),dat))))
    expect_equal(sum(fitted(f,type='full')['unused',,,,drop=FALSE]),0,tolerance=1e-12)
  }
  g <- ll_graph('path'); dat <- ll_data(g); dat <- dat[dat$B=='hi',,drop=FALSE]
  f <- fit_loglinear_graph(dat,g)
  full <- fitted(f,type='full')
  expect_false(anyNA(full))
  expect_equal(sum(full),1,tolerance=1e-12)
  expect_equal(sum(full[,'lo',,]),0)
})

test_that('Generators are deterministic maximal local cliques', {
  g <- ll_graph('mixed'); generators <- .loglinear_generators(g)
  expect_identical(generators,.loglinear_generators(g))
  expect_true(all(vapply(generators,function(s) .is_complete(g,s),logical(1))))
  expect_false(any(vapply(seq_along(generators),function(i) any(vapply(
    generators[-i],function(s) all(generators[[i]]%in%s),logical(1))),logical(1))))
})

test_that('Large state spaces remain factorized and full materialization is guarded', {
  p <- 70L; g <- igraph::set_vertex_attr(igraph::make_ring(p,circular=FALSE),'name',value=paste0('X',seq_len(p)))
  dat <- as.data.frame(stats::setNames(rep(list(factor(rep(c('a','b'),10))),p),igraph::V(g)$name))
  f <- fit_loglinear_graph(dat,g,engine='loglin')
  expect_equal(f$cells$full,2^70)
  expect_equal(f$cells$total_local,4*(p-1))
  expect_equal(f$cells$max_local,4)
  expect_true(is.finite(as.numeric(logLik(f))))
  expect_length(fitted(f)$atoms,p-1)
  expect_error(fitted(f,type='full'),'supported array size')
  expect_equal(sum(f$suffstats$observed_counts),20)
})

test_that('Full output guard is explicit and cannot override representational overflow', {
  g <- ll_graph('path'); f <- fit_loglinear_graph(ll_data(g),g,control=list(full_table_max_cells=4))
  expect_error(fitted(f,type='full'),'full_table_max_cells')
  expect_equal(sum(fitted(f,type='full',force=TRUE)),1,tolerance=1e-6)
  f$cells$full <- Inf
  expect_error(fitted(f,type='full',force=TRUE),'supported array size')
})

test_that('gRim and MASS independently agree with the loglin workflow', {
  skip_if_not_installed('gRim'); skip_if_not_installed('MASS')
  g <- ll_graph('mixed'); dat <- ll_data(g)
  a <- fit_loglinear_graph(dat,g,engine='loglin'); b <- fit_loglinear_graph(dat,g,engine='gRim')
  expect_equal(logLik(a),logLik(b),tolerance=1e-6)
  expect_equal(deviance(a),deviance(b),tolerance=1e-6)
  for (i in seq_along(a$local_fits)) expect_equal(a$local_fits[[i]]$fitted_prob,b$local_fits[[i]]$fitted_prob,tolerance=1e-6)
  oracle <- MASS::loglm(~A*B+B*C+C*D+D*A+C*D*E+E*F,data=do.call(table,dat),fit=TRUE,eps=1e-10,iter=10000)
  expect_equal(deviance(a),as.numeric(oracle$lrt),tolerance=1e-6)
  expect_equal(as.vector(fitted(a,type='full')),as.vector(oracle$fit/sum(oracle$fit)),tolerance=1e-6)
})

test_that('Invalid engine results and nonconvergence prevent successful composition', {
  g <- ll_graph('prime'); dat <- ll_data(g)
  local_mocked_bindings(.fit_loglinear_loglin=function(observed,...) list(fitted_counts=observed*0,
    iterations=NA_integer_,message='bad fit'))
  e <- tryCatch(fit_loglinear_graph(dat,g,engine='loglin'),atom_loglinear_error=identity)
  expect_s3_class(e,'atom_loglinear_error'); expect_identical(e$status,'invalid_counts')
  expect_equal(e$atom_id,1); expect_identical(e$engine,'loglin')
})

test_that('Generator violations and iteration exhaustion are distinct failures', {
  g <- ll_graph('prime'); dat <- ll_data(g)
  e <- tryCatch(fit_loglinear_graph(dat,g,engine='loglin',
    control=list(maxit=1L,tol=1e-15,contract_tol=1e-12)),atom_loglinear_error=identity)
  expect_s3_class(e,'atom_loglinear_error')
  expect_identical(e$status,'invalid_generator_margins')
  local_mocked_bindings(.fit_loglinear_loglin=function(observed,...) list(
    fitted_counts=array(sum(observed)/length(observed),dim=dim(observed),dimnames=dimnames(observed)),
    iterations=NA_integer_,message=''))
  e <- tryCatch(fit_loglinear_graph(dat,g,engine='loglin'),atom_loglinear_error=identity)
  expect_s3_class(e,'atom_loglinear_error')
  expect_identical(e$status,'invalid_generator_margins')
})

test_that('Explicit engine convergence warnings retain nonconvergence status', {
  local_mocked_bindings(loglin=function(...) { warning('algorithm did not converge'); list() },.package='stats')
  g <- ll_graph('prime')
  e <- tryCatch(fit_loglinear_graph(ll_data(g),g,engine='loglin'),atom_loglinear_error=identity)
  expect_s3_class(e,'atom_loglinear_error')
  expect_identical(e$status,'nonconverged')
})

test_that('Missing factors and duplicate declared levels are rejected', {
  g <- ll_graph('path'); dat <- ll_data(g)
  bad <- dat; bad$A[1] <- NA
  expect_error(fit_loglinear_graph(bad,g),'missing')
  bad <- dat; bad$A <- structure(as.integer(bad$A),levels=c('a','a'),class='factor')
  expect_error(fit_loglinear_graph(bad,g),'unique')
})

test_that('The default engine is gRim, engine choice is strict, and decomposition timing is isolated', {
  g <- ll_graph('path'); dat <- ll_data(g); dec <- atom_decompose(g)

  # default is equivalent to engine = 'gRim'
  a <- fit_loglinear_graph(dat, g, decomp = dec)
  b <- fit_loglinear_graph(dat, g, decomp = dec, engine = 'gRim')
  expect_identical(as.numeric(logLik(a)), as.numeric(logLik(b)))
  expect_identical(deviance(a), deviance(b))
  expect_identical(fitted(a), fitted(b))
  expect_identical(a$engine, 'gRim')
  expect_true(all(vapply(a$local_fits, function(z) z$requested_engine == 'gRim', logical(1))))

  # reusing a supplied decomposition costs zero decomposition time
  expect_identical(a$timing$decomposition, 0)

  # with decomp = NULL only atom_decompose() is timed
  c0 <- fit_loglinear_graph(dat, g)
  expect_true(is.numeric(c0$timing$decomposition))
  expect_length(c0$timing$decomposition, 1L)
  expect_true(is.finite(c0$timing$decomposition))
  expect_true(c0$timing$decomposition >= 0)
  expect_identical(as.numeric(logLik(c0)), as.numeric(logLik(a)))
  expect_identical(c0$decomp$atoms_index, dec$atoms_index)

  # an explicit engine never falls back: a failing gRim must surface as an error
  # on a graph whose atom is not complete, so the engine is really required
  gp <- ll_graph('prime'); datp <- ll_data(gp)
  local_mocked_bindings(.fit_loglinear_grim = function(...) stop('gRim unavailable'))
  e <- tryCatch(fit_loglinear_graph(datp, gp), atom_loglinear_error = identity)
  expect_s3_class(e, 'atom_loglinear_error')
  expect_identical(e$engine, 'gRim')
  expect_identical(e$attempted_engines, 'gRim')
  expect_match(conditionMessage(e), 'gRim unavailable', fixed = TRUE)

  # the removed auto option is rejected by match.arg
  expect_error(fit_loglinear_graph(dat, g, engine = 'auto'), 'arg')
})
