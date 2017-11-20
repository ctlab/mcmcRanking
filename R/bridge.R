mcmc_subgraph <- function(graph, module_size, iter) {
  edges <- data.frame(as_edgelist(graph, names = F)[,1:2] - 1)
  colnames(edges) <- c("from", "to")
  args <- c(nodes_size=length(V(graph)), module_size=module_size, iter=iter)
  res <- mcmc_subgraph_internal(edges, args) + 1
  return(V(graph)$name[res])
}

mcmc_sample <- function(graph, module_size, iter, times = 1) {
  edges <- data.frame(as_edgelist(graph, names = F)[,1:2] - 1)
  colnames(edges) <- c("from", "to")
  nodes <- data.frame(name=as.vector(V(graph)) - 1, likelihood = V(graph)$likelihood)
  args <- c(module_size=module_size, iter=iter, times=times)
  res <- mcmc_sample_internal(edges, nodes, args) + 1
  mat <- matrix(V(graph)$name[res], ncol = module_size, byrow = T)
  return(mat)
}

mcmc_onelong <- function(graph, module_size, start = 1, end) {
  edges <- data.frame(as_edgelist(graph, names = F)[,1:2] - 1)
  colnames(edges) <- c("from", "to")
  nodes <- data.frame(name=as.vector(V(graph)) - 1, likelihood = V(graph)$likelihood)
  args <- c(module_size=module_size, start=start, end=end)
  res <- mcmc_onelong_internal(edges, nodes, args)+1
  mat <- matrix(V(graph)$name[res], ncol = module_size, byrow = T)
  return(mat)
}

get_frequency <- function(graph, mat, inds = seq_len(nrow(mat))){
  prob <- table(V(graph)[mat[inds,]]$name)
  x <- numeric(length(V(graph)) - length(prob))
  names(x) <- setdiff(V(graph)$name, names(prob))
  return(c(prob, x))
}

get_prob <- function(graph, mat, inds = seq_len(nrow(mat))){
  return(get_frequency(graph, mat, inds)/length(inds))
}

get_inverse_likelihood <- function(graph, mat, inds = seq_len(nrow(mat))){
  inverse_likelihood <- 1 / apply(matrix(V(graph)[mat[inds,]]$likelihood, ncol = ncol(mat), byrow = T), 1, prod)
  inverse_likelihood <- sort(inverse_likelihood)
  return(sum(inverse_likelihood)/length(inds))
}
