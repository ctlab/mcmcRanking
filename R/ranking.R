#' Probabilistic ranking.
#'
#' Heuristic ranking maximizing the area under the ROC curve.
#'
#' @param graph An object of type \code{igraph}.
#' @param q Vector of probabilities of vertices to not belong to active module.
#' @return ranking.
#' @seealso \code{\link{mcmc_sample}}
#' @import igraph
#' @export
#' @examples
#' library(igraph)
#' data(exampleGraph)
#' rank <- probabilistic_rank(exampleGraph, V(exampleGraph)$q)
#' head(rank)
probabilistic_rank <- function(graph, q) {
  edgelist <- as_edgelist(graph, names = FALSE) - 1
  nodes <-
    data.frame(name = as.vector(V(graph)) - 1, q = q[V(graph)$name])
  res <- probabilistic_rank_internal(edgelist, nodes)
  names(res) <- V(graph)$name
  return(res)
}
