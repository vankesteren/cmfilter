context("[R] backend Product of Coefficients")
test_that("Single-core prodcoef cmf works", {
  o <- capture_output(
    res <<- cmf(
      d, 
      nStarts = 25,
      decisionFunction = cmfilter:::prodCoef,
      nCores = 1
    )
  )
  expect(inherits(res, "cmf"), "Result is not of class CMF")
})
test_that("Multi-core prodcoef cmf works", {
  o <- capture_output(
    res <<- cmf(
      d, 
      nStarts = 100,
      decisionFunction = cmfilter:::prodCoef,
      nCores = parallel::detectCores()
    )
  )
  expect(inherits(res, "cmf"), "Result is not of class CMF")
})
test_that("Update method works", {
  oldNstarts <- res$call$nStarts
  o <- capture_output(res <- update(res, 25))
  expect_equal(res$call$nStarts, oldNstarts + 25)
})
test_that("Print and summary methods work", {
  ptest <- capture_output_lines(print(res))
  expect_equal(ptest[2], "CMF Algorithm Results")
  stest <- capture_output_lines(summary(res))
  expect_equal(stest[2], "CMF Algorithm Results")
})
test_that("Screeplot method works", {
  fn <- tempfile(fileext = ".png")
  png(fn)
  screeplot(res, topn = 100)
  dev.off()
  expect_gt(file.size(fn), 318)
})
test_that("Plot method works", {
  fn <- tempfile(fileext = ".png")
  png(fn)
  plot(res)
  dev.off()
  expect_gt(file.size(fn), 318)
})
test_that("Cutoff setting works", {
  res <- setCutoff(res, cutoff = 0.1)
  expect_equal(res$call$cutoff, .1)
})
test_that("Adding method works", {
  res1 <- res
  res2 <- res
  res3 <- res1 + res2
  expect_equal(res3$selectionRate, (res1$selectionRate + res2$selectionRate)/2)
})


context("[R] backend Partial Correlation")
test_that("Single-core partcor cmf works", {
  o <- capture_output(
    res <<- cmf(
      d, 
      nStarts = 25,
      decisionFunction = cmfilter:::corMinusPartCor,
      nCores = 1
    )
  )
  expect(inherits(res, "cmf"), "Result is not of class CMF")
})
test_that("Multi-core partcor cmf works", {
  o <- capture_output(
    res <<- cmf(
      d, 
      nStarts = 100,
      decisionFunction = cmfilter:::corMinusPartCor,
      nCores = parallel::detectCores()
    )
  )
  expect(inherits(res, "cmf"), "Result is not of class CMF")
})
test_that("Update method works", {
  oldNstarts <- res$call$nStarts
  o <- capture_output(res <- update(res, 25))
  expect_equal(res$call$nStarts, oldNstarts + 25)
})
test_that("Print and summary methods work", {
  ptest <- capture_output_lines(print(res))
  expect_equal(ptest[2], "CMF Algorithm Results")
  stest <- capture_output_lines(summary(res))
  expect_equal(stest[2], "CMF Algorithm Results")
})
test_that("Screeplot method works", {
  fn <- tempfile(fileext = ".png")
  png(fn)
  screeplot(res, topn = 100)
  dev.off()
  expect_gt(file.size(fn), 318)
})
test_that("Plot method works", {
  fn <- tempfile(fileext = ".png")
  png(fn)
  plot(res)
  dev.off()
  expect_gt(file.size(fn), 318)
})
test_that("Cutoff setting works", {
  res <- setCutoff(res, cutoff = 0.1)
  expect_equal(res$call$cutoff, .1)
})
test_that("Adding method works", {
  res1 <- res
  res2 <- res
  res3 <- res1 + res2
  expect_equal(res3$selectionRate, (res1$selectionRate + res2$selectionRate)/2)
})


context("[R] backend Custom Selection Function")
test_that("Single-core custom cmf works", {
  # some arbitrary selection function slightly related to mediation
  selFun <<- function(x, m, y, crit = 3.84) {
    n   <- length(x)
    cxm <- crossprod(x, m) / (n - 1)
    cmy <- crossprod(m, y) / (n - 1)
    cxy <- crossprod(x, y) / (n - 1)
    tot <- cxm + cmy + cxy
    q <- abs(cxm * cmy / tot)
    return(q > crit)
  }
  
  o <- capture_output(
    res <<- cmf(
      d, 
      nStarts = 25,
      decisionFunction = selFun,
      nCores = 1,
      crit = 3.84
    )
  )
  expect(inherits(res, "cmf"), "Result is not of class CMF")
})
test_that("Multi-core custom cmf works", {
  o <- capture_output(
    res <<- cmf(
      d, 
      nStarts = 100,
      decisionFunction = selFun,
      nCores = parallel::detectCores(),
      crit = 3.84
    )
  )
  expect(inherits(res, "cmf"), "Result is not of class CMF")
})
test_that("Update method works", {
  oldNstarts <- res$call$nStarts
  o <- capture_output(res <- update(res, 25))
  expect_equal(res$call$nStarts, oldNstarts + 25)
})
test_that("Print and summary methods work", {
  ptest <- capture_output_lines(print(res))
  expect_equal(ptest[2], "CMF Algorithm Results")
  stest <- capture_output_lines(summary(res))
  expect_equal(stest[2], "CMF Algorithm Results")
})
test_that("Screeplot method works", {
  fn <- tempfile(fileext = ".png")
  png(fn)
  screeplot(res, topn = 100)
  dev.off()
  expect_gt(file.size(fn), 318)
})
test_that("Plot method works", {
  fn <- tempfile(fileext = ".png")
  png(fn)
  plot(res)
  dev.off()
  expect_gt(file.size(fn), 318)
})
test_that("Cutoff setting works", {
  res <- setCutoff(res, cutoff = 0.1)
  expect_equal(res$call$cutoff, .1)
})
test_that("Adding method works", {
  res1 <- res
  res2 <- res
  res3 <- res1 + res2
  expect_equal(res3$selectionRate, (res1$selectionRate + res2$selectionRate)/2)
})
