#' @export
mcmc <- function(mat, name) {
  if (!is.matrix(mat))
    stop("mat must be matrix.")
  if (!is.logical(mat))
    stop("mat mus be boolean vector.")
  if (!is.character(name))
    stop("name must be character vector.")
  structure(list(mat = mat,
                 name = name), class = "MCMC")
}



check_arguments <- function(graph, module_size, niter) {
  if (module_size > gorder(graph) || module_size < 1)
    stop("Required module size must be positive and not greather than graph size.")
  if (niter < 0)
    stop("number of iteration must be positive a number.")
  if (!("name" %in% vertex_attr_names(graph)))
    stop("graph must have 'name' attribute on vertices.")
}



#' Connected subgraph from uniform distribution.
#'
#' Generates a connected subgraph using Markov chain Monte Carlo (MCMC) method.
#'
#' @param graph An object of type \code{igraph} with \code{lieklihood} field in vertices.
#' @param module_size The size of subgraph.
#' @param niter Number of iterations.
#' @return Vector of vertex names of connected subgraph.
#' @seealso \code{\link{mcmc_sample}, \link{mcmc_onelong}}
#' @import igraph
#' @export
sample_subgraph <- function(graph, module_size, niter) {
  check_arguments(graph, module_size, niter)
  if (module_size == 0)
    return(c())
  edgelist <- as_edgelist(graph, names = F) - 1
  res <-
    sample_subgraph_internal(edgelist, gorder(graph), module_size, niter)
  return(V(graph)$name[which(res)])
}



#' @export
sample_llh <-
  function(graph,
           module_size = NULL,
           start_module = NULL,
           niter,
           exp_lh = 1,
           fixed_size = FALSE) {
    if (!xor(is.null(module_size) &&
             is.null(times),
             is.null(start_module))) {
      stop("One of the arguments module_size and times or start_module must be set.")
    }

    check_arguments(graph, module_size, niter)
    edgelist <- as_edgelist(graph, names = F) - 1

    if (!is.null(start_module)) {
      module_size <- sum(start_module[1, ])
    } else{
      start_module <-
        t(sample_subgraph_internal(edgelist, gorder(graph), module_size, 1))
    }

    res1 <-
      sample_llh_internal(edgelist,
                          V(graph)$likelihood ^ exp_lh,
                          niter,
                          fixed_size,
                          start_module)
    names(res1) <- seq_len(niter)
    return(res1)
  }



#' Set of connected subgraphs from uniform distribution.
#'
#' Generates set of independent subgraphs using Markov chain Monte Carlo (MCMC) method.
#'
#' @inheritParams sample_subgraph
#' @param times Number of subgraphs.
#' @param previous_mcmc Object of type MCMC.
#' @return Object of type MCMC.
#' @seealso \code{\link{sample_subgraph}, \link{mcmc_onelong}}
#' @import igraph
#' @export
#' @examples
#' graph <- barabasi.game(20, 1.2, 2, directed = F)
#' V(graph)$name <- letters[1:20]
#' V(graph)$likelihood <- 1
#' x <- mcmc_sample(graph = graph, module_size = 5, iter = 10, times = 4)
#' matrix(x$name[x$mat], nrow(x$mat))
mcmc_sample <-
  function(graph,
           module_size = NULL,
           times = NULL,
           previous_mcmc = NULL,
           niter,
           exp_lh,
           fixed_size = FALSE) {
    if (!xor(is.null(module_size) &&
             is.null(times),
             is.null(previous_mcmc))) {
      stop("One of the arguments module_size and times or previous_mcmc must be set.")
    }
    if (!is.null(previous_mcmc)) {
      if (class(previous_mcmc) != "MCMC")
        stop("previous_mcmc must be type of \"MCMC\"")
      start_module <- previous_mcmc$mat
    }
    check_arguments(graph, module_size, niter)
    edgelist <- as_edgelist(graph, names = F) - 1

    if (!is.null(previous_mcmc)) {
      module_size <- sum(start_module[1,])
      times <- nrow(start_module)
    } else{
      start_module <- t(replicate(
        times,
        sample_subgraph_internal(edgelist, gorder(graph), module_size, 1)
      ))
    }

    for (i in seq_along(exp_lh)) {
      res1 <- mcmc_sample_internal(edgelist,
                                   V(graph)$likelihood ^ exp_lh[i],
                                   fixed_size,
                                   niter,
                                   start_module)
      start_module <-  matrix(res1, ncol = gorder(graph), byrow = T)
    }

    return(mcmc(start_module, V(graph)$name))
  }



#' Set of connected subgraphs.
#'
#' Generates set of subgraphs using Markov chain Monte Carlo (MCMC) method.
#'
#' @param graph An object of type \code{igraph} with \code{lieklihood} field in vertices.
#' @param module_size The size of subgraph.
#' @param start Starting with this iteration, we write down all states of the Markov process.
#' @param end Number of iterations.
#' @return Object of type MCMC.
#' @seealso \code{\link{sample_subgraph}, \link{mcmc_sample}}
#' @import igraph
#' @export
mcmc_onelong <- function(graph, module_size, start, end) {
  check_arguments(graph, module_size, end)
  edgelist <- as_edgelist(graph, names = F) - 1
  res <-
    mcmc_onelong_internal(edgelist, V(graph)$likelihood, module_size, start, end) + 1
  ret <-
    mcmc(matrix(res, ncol = module_size, byrow = T), V(graph)$name)
  return(ret)
}



#' Frequency of vertices.
#'
#' The frequency of occurrence of vertices in one long run of method Markov chain Monte Carlo (MCMC).
#'
#' @param graph An object of type \code{igraph} with \code{lieklihood} field in vertices.
#' @param module_size The size of subgraph.
#' @param start Starting with this iteration, we write down all states of the Markov process.
#' @param end Number of iterations.
#' @return Named frequency vector.
#' @seealso \code{\link{mcmc_onelong}, \link{sample_subgraph}, \link{mcmc_sample}}
#' @import igraph
#' @export
mcmc_onelong_frequency <- function(graph, module_size, start, end) {
  check_arguments(graph, module_size, end)
  edgelist <- as_edgelist(graph, names = F) - 1
  res <-
    mcmc_onelong_frequency_internal(edgelist, V(graph)$likelihood, module_size, start, end)
  names(res) <- V(graph)$name
  return(res)
}
