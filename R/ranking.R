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
probabilistic_rank <- function(graph, q) {
  edgelist <- as_edgelist(graph, names = F) - 1
  nodes <-
    data.frame(name = as.vector(V(graph)) - 1, q = q[V(graph)$name])
  res <- probabilistic_rank_internal(edgelist, nodes)
  names(res) <- V(graph)$name
  return(res)
}



#' Ranking vertices.
#'
#' Ranking vertices according to the score.
#'
#' @param graph An object of type \code{igraph}.
#' @param q score.
#' @param size The initial size.
#' @return Named rank vector.
#' @seealso \code{\link{get_reverse_rank}}
#' @export
get_rank <- function(graph, q, size) {
  ret <- integer(gorder(graph))
  names(ret) <- V(graph)$name

  sorted_q <- names(sort(q))
  comps <-
    components(induced_subgraph(graph, V(graph)[sorted_q[seq_len(size)]]))
  best_comp_id <- order(comps$csize, decreasing = T)[1]
  bc <- names(which(comps$membership == best_comp_id))
  best <- bc[1]
  ret[bc] <- 1
  i <- 1
  while (size + i <= gorder(graph)) {
    comps <-
      components(induced_subgraph(graph, V(graph)[sorted_q[seq_len(size + i)]]))
    bc_new <-
      names(which(comps$membership == comps$membership[best]))
    i <- i + 1
    ret[setdiff(bc_new, bc)] <- max(ret) + 1
    bc <- bc_new
  }
  ret
}



#' Ranking vertices.
#'
#' Ranking vertices according to the score.
#'
#' @param graph An object of type \code{igraph}.
#' @param q score.
#' @return Named rank vector.
#' @seealso \code{\link{get_rank}}
#' @export
get_reverse_rank <- function(graph, q) {
  ret <- integer(gorder(graph))
  names(ret) <- V(graph)$name

  sorted_q <- names(sort(q, decreasing = T))
  cur_id <- seq_along(sorted_q)
  i <- 1
  while (length(cur_id) != 1) {
    comps <-
      components(induced_subgraph(graph, V(graph)[sorted_q[cur_id[-1]]]))
    best_comp_id <- order(comps$csize, decreasing = T)[1]
    bc <- names(which(comps$membership == best_comp_id))
    cur_new_id <- which(sorted_q %in% bc)
    ret[sorted_q[setdiff(cur_id, cur_new_id)]] <- i
    cur_id <- sort(cur_new_id)
    i <- i + 1
  }
  ret[sorted_q[cur_id]] <- i
  i + 1 - ret
}
