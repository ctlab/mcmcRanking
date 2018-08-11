#' Number of runs of iterative likelihood increase of vertices.
#'
#' Running MCMC with likelihood weights on vertices can cause to get stuck in locally
#' good vertices, thereby finding locally significant solutions.
#' So, we will iteratively increase vertex weights and this function calculates number of necessary iterations.
#'
#' @param x According to this value repetition depth calculates.
#' @return A number of runs.
#' @export
repetition_depth <- function(x) {
  d <- 0
  while (x > 4) {
    x <- sqrt(x)
    d <- d + 1
  }
  return(d)
}



#' Frequency of vertices.
#'
#' Calculates the frequency of occurrences of vertices in matrix object.
#'
#' @param mcmcObj Object of class MCMC.
#' @param inds Index numbers of rows involved in the calculation.
#' @return A named vector of frequency.
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
#' Set \emph{likelihood} vertex attribute to an \code{igraph} graph using pval vertex attribute.
#'
#' @param graph An \code{igraph} graph.
#' @param fdr Numeric constant, from the false discovery rate a p-value threshold is calculated.
#' @return An \code{igraph} graph with \emph{likelihood} vertex attribute.
#' @import igraph
#' @importFrom BioNet fitBumModel scoreFunction
#' @export
set_likelihood <- function(graph, fdr) {
  pvals <- V(graph)$pval
  names(pvals) <- V(graph)$name
  fb <- fitBumModel(pvals, plot = FALSE)
  V(graph)$likelihood <- exp(scoreFunction(fb = fb, fdr = fdr))
  graph
}
