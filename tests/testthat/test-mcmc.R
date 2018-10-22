context("mcmc")

test_that("bad call", {
  expect_error(mcmc(c(TRUE, TRUE, FALSE, TRUE)))
  expect_error(mcmc(matrix(c(1, 1, 0, 1), nrow = 1)))
  expect_error(mcmc(matrix(c(TRUE, TRUE, FALSE, TRUE),nrow = 1)))
})

test_that("mcmc constructor works", {
  mat <- matrix(c(TRUE, TRUE, FALSE, TRUE), 1, dimnames = list(c(), letters[1:4]))
  x <- mcmc(mat = mat)
  expect_identical(x$mat, mat)
})
