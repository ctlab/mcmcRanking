#include <Rcpp.h>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include "mcmc.h"
using namespace Rcpp;
using namespace std;
using mcmc::Graph;

vector<vector<unsigned>> make_edges(IntegerVector from, IntegerVector to, unsigned size){
  vector<vector<unsigned>> edges(size);
  for(int i = 0; i < from.size(); ++i){
    edges[from[i]].push_back(to[i]);
    edges[to[i]].push_back(from[i]);
  }
  return edges;
}

// [[Rcpp::export]]
IntegerVector mcmc_subgraph_internal(DataFrame df_edges, IntegerVector args) {
  IntegerVector from = df_edges["from"];
  IntegerVector to   = df_edges["to"];
  vector<double> nodes(args["nodes_size"]);
  vector<vector<unsigned>> edges = make_edges(from, to, args["nodes_size"]);
  Graph g = Graph(nodes, edges);
  g.initialize_module(g.random_subgraph(args["module_size"]));
  g.subgraph_iteration(args["iter"]);
  vector<unsigned> ret = g.get_inner_nodes();
  IntegerVector ret_(ret.begin(), ret.end());
  return ret_;
}

// [[Rcpp::export]]
IntegerVector mcmc_sample_internal(DataFrame df_edges, DataFrame df_nodes, IntegerVector args) {
  IntegerVector from = df_edges["from"];
  IntegerVector to   = df_edges["to"];
  IntegerVector names = df_nodes["name"];
  NumericVector likelihood = df_nodes["likelihood"];
  vector<double> nodes(likelihood.begin(), likelihood.end());
  vector<vector<unsigned>> edges = make_edges(from, to, names.size());
  Graph g = Graph(nodes, edges);
  vector<unsigned> module = g.random_subgraph(args["module_size"]);
  vector<unsigned> ret = g.sample_iteration(module, args["times"], args["iter"]);
  IntegerVector ret_(ret.begin(), ret.end());
  return ret_;
}

// [[Rcpp::export]]
IntegerVector mcmc_onelong_internal(DataFrame df_edges, DataFrame df_nodes, IntegerVector args) {
  IntegerVector from = df_edges["from"];
  IntegerVector to   = df_edges["to"];
  IntegerVector names = df_nodes["name"];
  NumericVector likelihood = df_nodes["likelihood"];
  vector<double> nodes(likelihood.begin(), likelihood.end());
  vector<vector<unsigned>> edges = make_edges(from, to, names.size());
  Graph g = Graph(nodes, edges);
  vector<unsigned> module = g.random_subgraph(args["module_size"]);
  g.initialize_module(module);
  vector<unsigned> ret = g.onelong_iteration(args["start"], args["end"]);
  IntegerVector ret_(ret.begin(), ret.end());
  return ret_;
}


// [[Rcpp::export]]
NumericVector mcmc_inverse_likelihood_internal(DataFrame df_edges, DataFrame df_nodes, IntegerVector args) {
  IntegerVector from = df_edges["from"];
  IntegerVector to   = df_edges["to"];
  NumericVector likelihood = df_nodes["likelihood"];
  vector<double> nodes(likelihood.begin(), likelihood.end());
  vector<vector<unsigned>> edges = make_edges(from, to, likelihood.size());
  Graph g = Graph(nodes, edges);
  vector<unsigned> module = g.random_subgraph(args["size"]);
  g.initialize_module(module);
  return(NumericVector::create(g.onelong_inverse_likelihood(args["start"], args["end"])));
}
