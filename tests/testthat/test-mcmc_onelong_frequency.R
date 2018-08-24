library(igraph)

context("mcmc_onelong_frequency")

test_that("bad call", {
  g <- make_full_graph(4)
  expect_error(mcmc_onelong_frequency(
    graph = g,
    module_size = 2,
    start = 50,
    niter = 100
  ))

  V(g)$name <- letters[1:4]
  V(g)$likelihood <- rep(1, 4)
  expect_error(mcmc_onelong_frequency(graph = g))
  expect_error(mcmc_onelong_frequency(graph = g, niter = 100))
  expect_error(mcmc_onelong_frequency(
    graph = g,
    module_size = 2,
    niter = 100
  ))
  expect_error(mcmc_onelong_frequency(
    graph = g,
    module_size = 2,
    start = 50
  ))
  expect_error(mcmc_onelong_frequency(
    graph = g,
    module_size = 5,
    start = 50,
    niter = 100
  ))
  expect_error(mcmc_onelong_frequency(
    graph = g,
    module_size = 2,
    start = 50,
    niter = -3
  ))
})

test_that("mcmc_sample works", {
  g <- make_full_graph(7)
  V(g)$name <- letters[1:7]
  V(g)$likelihood <- runif(7, 0.5, 2)
  freq <- mcmc_onelong_frequency(
    graph = g,
    module_size = 3,
    start = 50,
    niter = 100,
    fixed_size = TRUE
  )
  expect_identical(sum(freq), 150L)

  g <- make_full_graph(4)
  V(g)$name <- letters[1:4]
  V(g)$likelihood <- 4 ^ (2:-1)
  freq <-
    mcmc_onelong_frequency(
      graph = g,
      module_size = 2,
      start = 1e4,
      niter = 2e4,
      fixed_size = TRUE
    )
  p <- freq / 1e4
  expect_gte(p["a"], 84 / 89.25 - 0.05)
  expect_gte(p["b"], 69 / 89.25 - 0.05)
})
