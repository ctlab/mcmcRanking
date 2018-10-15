#' Probabilistic ranking.
#'
#' Heuristic ranking maximizing the area under the ROC curve.
#'
#' @param graph An object of type \code{igraph}.
#' @param q Vector of probabilities of vertices to not belong to active module.
#' @return ranking.
#' @seealso \code{\link{mcmc_sample}}
#' @importFrom  igraph as_edgelist V
#' @export
#' @examples
#' library(igraph)
#' data(exampleGraph)
#' q <- V(exampleGraph)$q
#' names(q) <- V(exampleGraph)$name
#' rank <- probabilistic_rank(exampleGraph, q)
#' head(rank)
probabilistic_rank <- function(graph, q) {
  if (is.null(names(q)))
    stop("Vector q must have a name attribute.")
  q <- q[V(graph)$name]
  if (any(is.na(q)))
    stop("Vector q does not contain probability for all vertices.")
  edgelist <- as_edgelist(graph, names = FALSE) - 1
  nodes <- data.frame(name = as.vector(V(graph)) - 1, q = q)
  res <- probabilistic_rank_internal(edgelist, nodes)
  names(res) <- V(graph)$name
  return(res)
}
