library(igraph)

context("sample_subgraph")

test_that("bad call", {
  g <- sample_gnm(10, 20)
  expect_error(sample_subgraph(
    graph = g,
    module_size = 5,
    niter = 1
  ))
})

test_that("sample_subgraph works", {
  g <- graph_from_atlas(100)
  V(g)$name <- letters[1:6]
  x <- sample_subgraph(g, 4, 1)
  expect_length(x, 4)
  expect_true(all(x %in% V(g)$name))

  g <- make_full_graph(20)
  V(g)$name <- letters[1:20]
  x <- sample_subgraph(g, 10, 1e2)
  expect_length(x, 10)
  expect_true(all(x %in% V(g)$name))

  # zero iteration is allowed
  x <- sample_subgraph(g, 5, 0)
  expect_length(x, 5)

  # returns empty vector if the module_size is zero
  x <- sample_subgraph(g, 0, 1)
  expect_length(x, 0)

  # returns all vertices from graph
  x <- sample_subgraph(g, 20, 10)
  expect_length(x, 20)

})
