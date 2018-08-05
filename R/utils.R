#' Number of runs of iterative likelihood increase of vertices.
#'
#' Running MCMC with likelihood weights on vertices can cause to get stuck in locally
#' good vertices, thereby finding locally significant solutions.
#' So, we will iteratively increase vertex weights and this function calculates number of necessary iterations.
#'
#' @param x According to this value repetition depth calculates.
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
#' @param mcmcObj Object of type MCMC.
#' @param inds Index numbers of rows involved in the calculation.
#' @return Named vector of frequency.
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
#' @param fdr Numeric constant, from the false discovery rate a p-value threshold is calculated.
#' @return igraph object.
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
