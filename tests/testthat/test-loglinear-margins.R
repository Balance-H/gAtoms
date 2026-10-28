test_that("Direct margins agree with small Cartesian-table oracles", {
  dat <- data.frame(A=factor(c('a','b','a','b'),levels=c('b','a','unused')),
    B=c('y','x','y','y'), C=c(2L,1L,2L,1L), D=c(TRUE,FALSE,FALSE,TRUE))
  sets <- list(c('A','B'),c('B','C'),c('A','B'),'D')
  z <- margin_counts(dat,sets)
  expect_s3_class(z,'loglinear_suffstats')
  expect_length(z$margins,3)
  expect_identical(z$levels$A,c('b','a','unused'))
  expect_identical(z$levels$C,c('1','2'))
  full <- do.call(table,Map(function(x,l) factor(x,levels=l),dat,z$levels))
  for (i in seq_along(z$sets)) {
    ref <- margin.table(full,match(z$sets[[i]],names(dat)))
    expect_equal(as.vector(z$margins[[i]]),as.vector(ref))
    expect_identical(dimnames(z$margins[[i]]),dimnames(ref))
  }
  expect_equal(z$cells$full,24)
  expect_equal(sum(z$observed_counts),nrow(dat))
  expect_equal(sum(z$margins[[1]]['unused',]),0)
})

test_that("Categorical input and local allocation limits are enforced", {
  dat <- data.frame(A=1:3,B=c('a','b','a'))
  expect_error(margin_counts(transform(dat,A=c(1.2,2,3)),list('A')),'integer-valued')
  for (bad in c(NA_real_,NaN,Inf,-Inf)) {
    x <- dat; x$A[1] <- bad
    expect_error(margin_counts(x,list('A')))
  }
  x <- dat; names(x) <- c('A','A')
  expect_error(margin_counts(x,list('A')),'unique')
  expect_error(margin_counts(dat,list(c('A','A'))),'distinct')
  expect_error(margin_counts(dat,list('unknown')),'known')
  expect_error(margin_counts(dat,list(c('A','B')),max_cells=2),'max_cells')
  expect_identical(margin_counts(data.frame(A=ordered(c('b','a'),levels=c('c','b','a'))),list('A'))$levels$A,c('c','b','a'))
  expect_identical(.loglinear_cells(rep(2,1100)),Inf)
})

test_that("Variable names cannot collide with configuration encoding arguments", {
  dat <- data.frame(sep=c('a','a','b'),collapse=c(1L,1L,2L))
  z <- margin_counts(dat,list(c('sep','collapse')))
  expect_identical(z$variables,c('sep','collapse'))
  expect_equal(sort(z$observed_counts),c(1L,2L))
})
