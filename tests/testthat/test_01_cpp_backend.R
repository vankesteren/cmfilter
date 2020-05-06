# Unit tests for cpp backend of CMF\


context("Data Generation")

# generate data
test_that("Data generation works", {
  set.seed(45)
  apaths <- runif(100, -1, 1)
  bpaths <- runif(100, -1, 1)
  d <<- generateMed(50, apaths, bpaths, r2y = .6, dir = F)
  expect_s3_class(d, "data.frame")
})
test_that("Complex data generation works", {
  set.seed(45)
  p <- 16
  P <- qr.Q(qr(matrix(rnorm(p^2), p))) # eigenvectors
  rate <- 1.1
  e <- (rate^(p:1)/rate*p)/sum(rate^(p:1)/rate) # eigenvalues sum to p
  S <- cov2cor(crossprod(P, P * e))
  apaths <- c(0.3, rep(0, 15))
  bpaths <- c(0.3, sign(S)[-1,1]*c(rep(0.8, 3), rep(0.4, 12)))
  Sigma <- diag(1 - apaths^2)
  S <- S * tcrossprod(diag(Sigma))
  diag(S) <- 0
  Sigma <- Sigma + S
  rsquared <- 0.5
  
  d <- generateMed(n = 400, 
                   a = apaths,
                   b = bpaths,
                   Sigma = Sigma,
                   residual = TRUE,
                   r2y = rsquared,
                   empirical = TRUE)
  
  expect_s3_class(d, "data.frame")
})


context("CPP backend Product of Coefficients")
test_that("Single-core prodcoef cmf works", {
  res <<- cmf(
    d, 
    nStarts = 100,
    decisionFunction = "prodcoef",
    nCores = 1,
    pb = FALSE
  )
  expect(inherits(res, "cmf"), "Result is not of class CMF")
})
test_that("Multi-core prodcoef cmf works", {
  res <<- cmf(
    d, 
    nStarts = 400,
    decisionFunction = "prodcoef",
    nCores = 2,
    pb = FALSE
  )
  expect(inherits(res, "cmf"), "Result is not of class CMF")
})
test_that("Update method works", {
  oldNstarts <- res$call$nStarts
  res <- update(res, 100)
  expect_equal(res$call$nStarts, oldNstarts + 100)
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
  expect_equal(res$cutoff, .1)
})
test_that("Adding method works", {
  res1 <- res
  res2 <- res
  res3 <- res1 + res2
  expect_equal(res3$selectionRate, (res1$selectionRate + res2$selectionRate)/2)
})


context("CPP backend Causal Steps")
test_that("Single-core csteps cmf works", {
  res <<- cmf(
    d, 
    nStarts = 100,
    decisionFunction = "causalsteps",
    nCores = 1,
    pb = FALSE
  )
  expect(inherits(res, "cmf"), "Result is not of class CMF")
})
test_that("Multi-core csteps cmf works", {
  res <<- cmf(
    d, 
    nStarts = 400,
    decisionFunction = "causalsteps",
    nCores = 2,
    pb = FALSE
  )
  expect(inherits(res, "cmf"), "Result is not of class CMF")
})
test_that("Update method works", {
  oldNstarts <- res$call$nStarts
  res <- update(res, 100)
  expect_equal(res$call$nStarts, oldNstarts + 100)
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
  expect_equal(res$cutoff, .1)
})
test_that("Adding method works", {
  res1 <- res
  res2 <- res
  res3 <- res1 + res2
  expect_equal(res3$selectionRate, (res1$selectionRate + res2$selectionRate)/2)
})