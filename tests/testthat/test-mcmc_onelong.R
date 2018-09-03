library(igraph)

context("mcmc_onelong")

test_that("bad call", {
  g <- make_full_graph(4)
  expect_error(mcmc_onelong(
    graph = g,
    subgraph_order = 2,
    start = 50,
    niter = 100
  ))

  V(g)$name <- letters[1:4]
  V(g)$likelihood <- rep(1, 4)
  expect_error(mcmc_onelong(graph = g))
  expect_error(mcmc_onelong(graph = g, niter = 100))
  expect_error(mcmc_onelong(
    graph = g,
    subgraph_order = 2,
    niter = 100
  ))
  expect_error(mcmc_onelong(
    graph = g,
    subgraph_order = 2,
    start = 50
  ))
  expect_error(mcmc_onelong(
    graph = g,
    subgraph_order = 5,
    start = 50,
    niter = 100
  ))
  expect_error(mcmc_onelong(
    graph = g,
    subgraph_order = 2,
    start = 50,
    niter = -3
  ))
})

test_that("mcmc_sample works", {
  g <- make_full_graph(7)
  V(g)$name <- letters[1:7]
  V(g)$likelihood <- runif(7, 0.5, 2)
  x <-
    mcmc_onelong(
      graph = g,
      subgraph_order = 3,
      start = 50,
      niter = 100,
      fixed_size = TRUE
    )
  expect_identical(sum(get_frequency(x)), 50 * 3)

  g <- make_full_graph(4)
  V(g)$name <- letters[1:4]
  V(g)$likelihood <- 4 ^ (2:-1)
  x <-
    mcmc_onelong(
      graph = g,
      subgraph_order = 2,
      start = 1e4,
      niter = 2e4,
      fixed_size = TRUE
    )
  p <- get_frequency(x) / 1e4
  expect_gte(p["a"], 84 / 89.25 - 0.05)
  expect_gte(p["b"], 69 / 89.25 - 0.05)
})
