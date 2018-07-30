library(igraph)

context("probabilistic_rank")

test_that("probabilistic_rank works", {
  g <- make_full_graph(20)
  V(g)$name <- letters[seq_len(20)]
  q <- 0.05 * seq_len(20)
  names(q) <- letters[seq_len(20)]
  x <- probabilistic_rank(g, q)
  expect_equivalent(x[letters[seq_len(20)]], seq_len(20))

  g <- graph_from_atlas(286)
  V(g)$name <- letters[seq_len(7)]
  q <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.1)
  names(q) <- V(g)$name
  x <- probabilistic_rank(g, q)
  expect_equivalent(x[letters[seq_len(7)]], c(2, 3, 4, 5, 6, 7, 1))

  q <- c(0.8, 0.8, 0.6, 0.2, 0.2, 0.4, 1)
  names(q) <- V(g)$name
  x <- probabilistic_rank(g, q)
  expect_equivalent(x[letters[seq_len(7)]], c(6, 5, 4, 2, 1, 3, 7))

  g <- graph_from_atlas(14)
  V(g)$name <- letters[seq_len(4)]
  q <- c(0.9, 0.5, 0.4, 0.1)
  names(q) <- letters[seq_len(4)]
  x <- probabilistic_rank(g, q)
  expect_equivalent(x[letters[seq_len(4)]], c(4, 2, 1, 4))
})
