#include <Rcpp.h>
#include <vector>
#include "mcmc.h"
#include "utils.h"
using namespace Rcpp;
using namespace std;
using mcmc::Graph;

vector<vector<unsigned>> adj_list(IntegerMatrix edgelist, size_t gorder){
  vector<vector<unsigned>> edges(gorder);
  for(int i = 0; i < edgelist.nrow(); ++i){
    edges[edgelist(i, 0)].push_back(edgelist(i, 1));
    edges[edgelist(i, 1)].push_back(edgelist(i, 0));
  }
  return edges;
}

// [[Rcpp::export]]
LogicalVector sample_subgraph_internal(IntegerMatrix edgelist, int gorder, int module_size, size_t niter) {
  vector<double> nodes(gorder, 1);
  Graph g = Graph(nodes, adj_list(edgelist, gorder), true);
  vector<vector<unsigned>> module;
  module.push_back(g.random_subgraph(module_size));
  vector<char> ret = g.sample_iteration(module, 1, niter);
  return LogicalVector (ret.begin(), ret.end());
}

// [[Rcpp::export]]
NumericVector sample_llh_internal(IntegerMatrix edgelist, NumericVector likelihood, size_t niter, bool fixed_size, LogicalMatrix start_module) {
  Graph g = Graph(likelihood, adj_list(edgelist, likelihood.size()), fixed_size);
  vector<unsigned> module;
  for(int j = 0; j < start_module.ncol(); ++j){
    if(start_module(0, j)){
      module.push_back(j);
    }
  }
  vector<double> ret = g.sample_llh(module, niter);
  return NumericVector (ret.begin(), ret.end());
}

// [[Rcpp::export]]
LogicalVector mcmc_sample_internal(IntegerMatrix edgelist, NumericVector likelihood, bool fixed_size, size_t niter, LogicalMatrix start_module) {
  Graph g = Graph(likelihood, adj_list(edgelist, likelihood.size()), fixed_size);
  vector<vector<unsigned>> module;
  for(int i = 0; i < start_module.nrow(); ++i){
    module.push_back(vector<unsigned>());
    for(int j = 0; j < start_module.ncol(); ++j){
      if(start_module(i, j)){
        module[i].push_back(j);
      }
    }
  }
  vector<char> ret = g.sample_iteration(module, module.size(), niter);
  return LogicalVector (ret.begin(), ret.end());
}

// [[Rcpp::export]]
IntegerVector mcmc_onelong_internal(IntegerMatrix edgelist, NumericVector likelihood, int module_size, size_t start, size_t end) {
  Graph g = Graph(likelihood, adj_list(edgelist, likelihood.size()), true);
  vector<unsigned> module = g.random_subgraph(module_size);
  g.initialize_module(module);
  vector<char> ret = g.onelong_iteration(start, end);
  return IntegerVector (ret.begin(), ret.end());
}

// [[Rcpp::export]]
IntegerVector mcmc_onelong_frequency_internal(IntegerMatrix edgelist, NumericVector likelihood, int module_size, size_t start, size_t end) {
  Graph g = Graph(likelihood, adj_list(edgelist, likelihood.size()), true);
  vector<unsigned> module = g.random_subgraph(module_size);
  g.initialize_module(module);
  vector<unsigned> ret = g.onelong_iteration_frequency(start, end);
  return IntegerVector (ret.begin(), ret.end());
}
