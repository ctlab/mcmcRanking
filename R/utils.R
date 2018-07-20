#' @export
repetition_depth <- function(x) {
  d <- 0
  while (x > 4) {
    x <- sqrt(x)
    d <- d + 1
  }
  return(d)
}



#' Frequency of vertecies.
#'
#' Calculates the frequency of occurences of vertices in matrix object.
#'
#' @param mcmcObj Object of type MCMC.
#' @param inds Index numbers of rows involved in the calculation.
#' @return Named vector of frequency.
#' @seealso \code{\link{get_prob}}
#' @import igraph
#' @export
get_frequency <-
  function(mcmcObj, inds = seq_len(nrow(mcmcObj$mat))) {
    freq <- colSums(mcmcObj$mat[inds, , drop = FALSE])
    names(freq) <- mcmcObj$name
    return(freq)
  }


#' Set likelihood.
#'
#' Set likelihood attribute to vertices of graph using pval attribute.
#'
#' @param graph An object of type \code{igraph}.
#' @return igraph object.
#' @import igraph
#' @import BioNet
#' @export
set_likelihood <- function(graph, fdr) {
  fb <- fitBumModel(V(graph)$pval, plot = FALSE)
  V(graph)$likelihood <- exp(scoreFunction(fb = fb, fdr = fdr))
  graph
}



#' Vertex probability.
#'
#' Accurate estimate of vertex probability using all connected subgraphs.
#'
#' @param graph An object of type \code{igraph} with \code{lieklihood} field in vertices.
#' @return Named vector of probabilites.
#' @details Time complexity of method is \strong{exponential}. Use it only for the graphs of size less than 30.
#' @seealso \code{\link{get_prob}}
#' @import igraph
#' @export
real_prob <- function(graph) {
  edgelist <- as_edgelist(graph, names = F) - 1
  res <- real_prob_internal(edgelist, V(graph)$likelihood)
  names(res) <- V(graph)$name
  return(res)
}
