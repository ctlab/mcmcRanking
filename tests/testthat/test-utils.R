context("utils")

test_that("repetition_depth works", {
  for (i in 0:2) {
    for (j in c(2.01, 3.1, 4)) {
      expect_equal(repetition_depth(j ^ (2 ^ i)), i)
    }
  }
})

test_that("get_frequency works", {
  mcmcobj <-
    mcmc(matrix(sample(c(TRUE, FALSE), 100, TRUE), 20), letters[1:5])
  ans <- colSums(mcmcobj$mat)
  names(ans) <- letters[1:5]
  freq <- get_frequency(mcmcobj)
  expect_identical(freq, ans)
  expect_setequal(names(freq), names(ans))
})

test_that("set_likelihood works", {
  data("exampleGraph")
  g <- igraph::delete_vertex_attr(exampleGraph, "likelihood")
  g <- set_likelihood(g, 1e-7)
  expect_true(any("likelihood" == igraph::vertex_attr_names(g)))
  expect_true(is.numeric(igraph::V(g)$likelihood))
})
