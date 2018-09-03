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
check_arguments <- function(graph, subgraph_order, niter) {
  if (subgraph_order > gorder(graph) || subgraph_order < 0)
    stop("Subgraph order must be non-negative and not greather than graph order")
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



#' Sample a connected subgraph from an uniform distribution.
#'
#' Generates a connected subgraph from an uniform distribution using Markov
#' chain Monte Carlo (MCMC) method.
#'
#' @param graph The original graph of class \code{igraph}.
#' @param subgraph_order The order of subgraph.
#' @param niter Number of iterations the MCMC method will run.
#' @return Character vector.
#' @seealso \code{\link{mcmc_sample}}
#' @importFrom igraph as_edgelist gorder V
#' @export
#' @examples
#' data(exampleGraph)
#' sample_subgraph(exampleGraph, 10, 1e4)
sample_subgraph <- function(graph, subgraph_order, niter) {
  check_arguments(graph, subgraph_order, niter)
  edgelist <- as_edgelist(graph, names = FALSE) - 1
  res <-
    sample_subgraph_internal(edgelist, gorder(graph), subgraph_order, niter)
  return(V(graph)$name[which(res)])
}



#' Sample log-likelihood values during one run.
#'
#' Generates log-likelihood values during one Markov chain Monte Carlo (MCMC)
#' run for analyzing behavior of subgraph.
#'
#' \code{graph} must contain vertex attribute \emph{likelihood}.
#'
#' @inheritParams sample_subgraph
#' @param exp_lh The power to which the likelihood values will be raised
#' @param fixed_order Logical scalar. Whether to fix the subgraph order.
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
           subgraph_order,
           niter,
           exp_lh = 1,
           fixed_order = FALSE) {
    if (missing(subgraph_order) && !fixed_order) {
      subgraph_order <- 0
    }
    check_arguments(graph, subgraph_order, niter)
    edgelist <- as_edgelist(graph, names = FALSE) - 1

    start_module <-
      t(sample_subgraph_internal(edgelist, gorder(graph), subgraph_order, 1))

    llhs <-
      sample_llh_internal(edgelist,
                          V(graph)$likelihood ^ exp_lh,
                          niter,
                          fixed_order,
                          start_module)
    return(setNames(llhs, seq_len(niter)))
  }



#' Sample connected subgraphs by likelihood value.
#'
#' Generates set of independent subgraphs using Markov chain Monte Carlo (MCMC)
#' method.
#'
#' @inherit sample_llh details params
#' @param times A number of sampling subgraphs.
#' @param previous_mcmc Object of class MCMC.
#' @return Object of class MCMC.
#' @seealso \code{\link{sample_llh}, \link{sample_subgraph},
#'   \link{mcmc_onelong}, \link{get_frequency}}
#' @importFrom igraph as_edgelist gorder V
#' @export
#' @examples
#' data(exampleGraph)
#' x <- mcmc_sample(exampleGraph, subgraph_order = 0, times = 1e3, niter = 100)
#' freq <- get_frequency(x)
#' tail(sort(freq))
mcmc_sample <-
  function(graph,
           subgraph_order,
           times,
           niter,
           previous_mcmc,
           exp_lh = 1,
           fixed_order = FALSE) {
    if (missing(previous_mcmc) && missing(subgraph_order) && !fixed_order) {
      subgraph_order <- 0
    }
    if (!xor(missing(subgraph_order) &&
             missing(times),
             missing(previous_mcmc))) {
      stop("subgraph_order and times or previous_mcmc must be specified.")
    }
    if (!missing(previous_mcmc)) {
      if (class(previous_mcmc) != "MCMC")
        stop("previous_mcmc must be class of \"MCMC\"")
      start_module <- previous_mcmc$mat
      subgraph_order <- sum(start_module[1,])
    }
    check_arguments(graph, subgraph_order, niter)
    edgelist <- as_edgelist(graph, names = FALSE) - 1

    if (missing(previous_mcmc)) {
      start_module <-
        matrix(unlist(
          replicate(
            times,
            sample_subgraph_internal(edgelist, gorder(graph), subgraph_order, 1),
            simplify = FALSE
          )
        ), nrow = times, byrow = TRUE)
    }
    likelihoods <- outer(V(graph)$likelihood, exp_lh, "^")
    res1 <- mcmc_sample_internal(edgelist,
                                 likelihoods,
                                 fixed_order,
                                 niter,
                                 start_module)
    start_module <- matrix(res1, ncol = gorder(graph), byrow = TRUE)
    return(mcmc(start_module, V(graph)$name))
  }



#' Sample connected subgraphs by likelihood value during one long MCMC run.
#'
#' Generates set of subgraphs using Markov chain Monte Carlo (MCMC) method.
#'
#' @inherit sample_llh details params
#' @param start Starting with this iteration, we write down all states of the
#'   Markov process.
#' @return Object of class MCMC.
#' @seealso \code{\link{mcmc_onelong_frequency}, \link{mcmc_sample},
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
           subgraph_order,
           start,
           niter,
           fixed_order = FALSE) {
    if (missing(subgraph_order) && !fixed_order) {
      subgraph_order <- 0
    }
    check_arguments(graph, subgraph_order, niter)
    edgelist <- as_edgelist(graph, names = FALSE) - 1
    res <-
      mcmc_onelong_internal(edgelist,
                            V(graph)$likelihood,
                            fixed_order,
                            subgraph_order,
                            start,
                            niter)
    ret <-
      mcmc(matrix(res, ncol = gorder(graph), byrow = TRUE), V(graph)$name)
    return(ret)
  }



#' The frequency of occurrence of vertices in one long MCMC run.
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
           subgraph_order,
           start,
           niter,
           fixed_order = FALSE) {
    if (missing(subgraph_order) && !fixed_order) {
      subgraph_order <- 0
    }
    check_arguments(graph, subgraph_order, niter)
    edgelist <- as_edgelist(graph, names = FALSE) - 1
    res <-
      mcmc_onelong_frequency_internal(edgelist,
                                      V(graph)$likelihood,
                                      fixed_order,
                                      subgraph_order,
                                      start,
                                      niter)
    names(res) <- V(graph)$name
    return(res)
  }
