library(igraph)

context("sample_llh")

test_that("sample_llh works", {
  g <- make_ring(5)
  V(g)$name <- letters[1:5]
  V(g)$likelihood <- 1:5
  x <-
    sample_llh(
      graph = g,
      subgraph_order = 2,
      niter = 100,
      fixed_size = TRUE
    )
  expect_length(x, 100)
  expect_true(all(x %in% log(c(2, 6, 12, 20, 5))))

  g <- make_full_graph(5)
  V(g)$name <- letters[1:5]
  V(g)$likelihood <- c(1e3, 1e3, 1, 1, 1)
  x <-
    sample_llh(
      graph = g,
      subgraph_order = 2,
      niter = 2e3,
      fixed_size = TRUE
    )
  expect_gte(sum(tail(x, 1e3) == log(1e6)), 950)

  # must return vector of zeros
  x <-
    sample_llh(
      graph = g,
      subgraph_order = 0,
      niter = 100,
      fixed_size = TRUE
    )
  expect_true(all(x == rep(0, 100)))
})
