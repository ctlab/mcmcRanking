#' MCMC class constructor
#'
#' MCMC data structure contains information about subgraph.
#'
#' @param mat A boolean matrix object where every row defines a subgraph.
#' @param name A character vector contains names of graph vertices.
#' @return Object of class MCMC.
#' @export
mcmc <- function(mat, name) {
  if (!is.matrix(mat))
    stop("mat must be a matrix.")
  if (!is.logical(mat))
    stop("mat must be a boolean vector.")
  if (!is.character(name))
    stop("name must be a character vector.")
  if (ncol(mat) != length(name))
    stop("Number of columns in mat is not equal to size of name.")
  structure(list(mat = mat,   name = name), class = "MCMC")
}


#' @importFrom igraph gorder vertex_attr_names is_simple
check_arguments <- function(graph, module_size, niter) {
  if (module_size > gorder(graph) || module_size < 0)
    stop("Module size must be non-negative and not greather than graph size.")
  if (niter < 0)
    stop("Number of iteration must be a non-negative number.")
  if (!("name" %in% vertex_attr_names(graph)))
    stop("graph must have \"name\" attribute on vertices.")
  if (!is_simple(graph))
    stop(
      "A simple graph is expected.\nSimple graphs are graphs which do not
      contain loop and multiple edges."
    )
}



#' Connected subgraph from uniform distribution.
#'
#' Generates a connected subgraph using Markov chain Monte Carlo (MCMC) method.
#'
#' @param graph An \code{igraph} graph with \code{lieklihood} vertex attribute.
#' @param module_size The size of subgraph.
#' @param niter Number of iterations.
#' @return Vector of vertex names of connected subgraph.
#' @seealso \code{\link{mcmc_sample}, \link{mcmc_onelong}}
#' @importFrom igraph as_edgelist gorder V
#' @export
#' @examples
#' data(exampleGraph)
#' sample_subgraph(exampleGraph, 10, 1e4)
sample_subgraph <- function(graph, module_size, niter) {
  check_arguments(graph, module_size, niter)
  edgelist <- as_edgelist(graph, names = FALSE) - 1
  res <-
    sample_subgraph_internal(edgelist, gorder(graph), module_size, niter)
  return(V(graph)$name[which(res)])
}



#' Generates log-likelihood values during one MCMC run.
#'
#' Generates log-likelihood values during one MCMC run for analyzing behavior of
#' subgraph.
#'
#' @inheritParams sample_subgraph
#' @param exp_lh The exponent likelihood values to be raised to.
#' @param fixed_size \code{TRUE} if the module size is fixed.
#' @return A named vector of likelihoods where names are number of iteration.
#' @importFrom igraph as_edgelist gorder V
#' @importFrom stats setNames
#' @export
#' @examples
#' data(exampleGraph)
#' llhs <- sample_llh(exampleGraph, 0, 1e5)
#' tail(llhs)
sample_llh <-
  function(graph,
           module_size,
           niter,
           exp_lh = 1,
           fixed_size = FALSE) {
    if (missing(module_size) && !fixed_size) {
      module_size <- 0
    }
    check_arguments(graph, module_size, niter)
    edgelist <- as_edgelist(graph, names = FALSE) - 1

    start_module <-
      t(sample_subgraph_internal(edgelist, gorder(graph), module_size, 1))

    llhs <-
      sample_llh_internal(edgelist,
                          V(graph)$likelihood ^ exp_lh,
                          niter,
                          fixed_size,
                          start_module)
    return(setNames(llhs, seq_len(niter)))
  }



#' Sampling set of connected subgraphs by likelihood value.
#'
#' Generates set of independent subgraphs using Markov chain Monte Carlo (MCMC)
#' method.
#'
#' @inheritParams sample_llh
#' @param times A number of subgraphs.
#' @param previous_mcmc Object of class MCMC.
#' @return Object of class MCMC.
#' @seealso \code{\link{sample_llh}, \link{sample_subgraph},
#'   \link{mcmc_onelong}, \link{get_frequency}}
#' @importFrom igraph as_edgelist gorder V
#' @export
#' @examples
#' data(exampleGraph)
#' x <- mcmc_sample(exampleGraph, module_size = 0, times = 1e3, niter = 100)
#' freq <- get_frequency(x)
#' tail(sort(freq))
mcmc_sample <-
  function(graph,
           module_size,
           times,
           niter,
           previous_mcmc,
           exp_lh = 1,
           fixed_size = FALSE) {
    if (missing(previous_mcmc) && missing(module_size) && !fixed_size) {
      module_size <- 0
    }
    if (!xor(missing(module_size) &&
             missing(times),
             missing(previous_mcmc))) {
      stop("module_size and times or previous_mcmc must be specified.")
    }
    if (!missing(previous_mcmc)) {
      if (class(previous_mcmc) != "MCMC")
        stop("previous_mcmc must be class of \"MCMC\"")
      start_module <- previous_mcmc$mat
      module_size <- sum(start_module[1,])
    }
    check_arguments(graph, module_size, niter)
    edgelist <- as_edgelist(graph, names = FALSE) - 1

    if (missing(previous_mcmc)) {
      start_module <-
        matrix(unlist(
          replicate(
            times,
            sample_subgraph_internal(edgelist, gorder(graph), module_size, 1),
            simplify = FALSE
          )
        ), nrow = times, byrow = TRUE)
    }
    for (i in seq_along(exp_lh)) {
      res1 <- mcmc_sample_internal(edgelist,
                                   V(graph)$likelihood ^ exp_lh[i],
                                   fixed_size,
                                   niter,
                                   start_module)
      start_module <-
        matrix(res1, ncol = gorder(graph), byrow = TRUE)
    }
    return(mcmc(start_module, V(graph)$name))
  }



#' Sampling set of connected subgraphs by likelihood value.
#'
#' Generates set of subgraphs using Markov chain Monte Carlo (MCMC) method.
#'
#' @param graph An \code{igraph} graph with \code{lieklihood} vertex attribute.
#' @param module_size The size of subgraph.
#' @param start Starting with this iteration, we write down all states of the
#'   Markov process.
#' @param niter Number of iterations.
#' @param fixed_size \code{TRUE} if the module size is fixed.
#' @return Object of class MCMC.
#' @seealso \code{\link{sample_subgraph}, \link{mcmc_sample},
#'   \link{get_frequency}}
#' @importFrom igraph as_edgelist gorder V
#' @export
#' @examples
#' data(exampleGraph)
#' x <- mcmc_onelong(exampleGraph, 50, 1e4, 2e4)
#' freq <- get_frequency(x)
#' tail(sort(freq))
mcmc_onelong <-
  function(graph,
           module_size,
           start,
           niter,
           fixed_size = FALSE) {
    if (missing(module_size) && !fixed_size) {
      module_size <- 0
    }
    check_arguments(graph, module_size, niter)
    edgelist <- as_edgelist(graph, names = FALSE) - 1
    res <-
      mcmc_onelong_internal(edgelist,
                            V(graph)$likelihood,
                            fixed_size,
                            module_size,
                            start,
                            niter)
    ret <-
      mcmc(matrix(res, ncol = gorder(graph), byrow = TRUE), V(graph)$name)
    return(ret)
  }



#' Frequency of vertices using one long run in MCMC.
#'
#' The frequency of occurrence of vertices in one long run of method Markov
#' chain Monte Carlo (MCMC).
#'
#' @inheritParams mcmc_onelong
#' @return A named frequency vector.
#' @seealso \code{\link{mcmc_onelong}, \link{sample_subgraph},
#'   \link{mcmc_sample}}
#' @importFrom igraph as_edgelist V
#' @export
#' @examples
#' data(exampleGraph)
#' freq <- mcmc_onelong_frequency(exampleGraph, 50, 1e4, 2e4)
#' tail(sort(freq), 60)
mcmc_onelong_frequency <-
  function(graph,
           module_size,
           start,
           niter,
           fixed_size = FALSE) {
    if (missing(module_size) && !fixed_size) {
      module_size <- 0
    }
    check_arguments(graph, module_size, niter)
    edgelist <- as_edgelist(graph, names = FALSE) - 1
    res <-
      mcmc_onelong_frequency_internal(edgelist,
                                      V(graph)$likelihood,
                                      fixed_size,
                                      module_size,
                                      start,
                                      niter)
    names(res) <- V(graph)$name
    return(res)
  }
