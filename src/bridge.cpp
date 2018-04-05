#include <Rcpp.h>
#include <vector>
#include <map>
#include <string>
#include <queue>
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
  vector<double> nodes(args["nodes_size"], 1);
  vector<vector<unsigned>> edges = make_edges(from, to, args["nodes_size"]);
  Graph g = Graph(nodes, edges);
  vector<unsigned> module = g.random_subgraph(args["module_size"]);
  vector<unsigned> ret = g.sample_iteration(module, 1, args["iter"]);
  IntegerVector ret_(ret.begin(), ret.end());
  return ret_;
}

// [[Rcpp::export]]
IntegerVector mcmc_sample_internal(DataFrame df_edges, DataFrame df_nodes, IntegerVector args, IntegerVector start_module) {
  IntegerVector from = df_edges["from"];
  IntegerVector to   = df_edges["to"];
  IntegerVector names = df_nodes["name"];
  NumericVector likelihood = df_nodes["likelihood"];
  vector<double> nodes(likelihood.begin(), likelihood.end());
  vector<vector<unsigned>> edges = make_edges(from, to, names.size());
  Graph g = Graph(nodes, edges);
  vector<unsigned> module(start_module.begin(), start_module.end());
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
IntegerVector mcmc_onelong_frequency_internal(DataFrame df_edges, DataFrame df_nodes, IntegerVector args) {
  IntegerVector from = df_edges["from"];
  IntegerVector to   = df_edges["to"];
  IntegerVector names = df_nodes["name"];
  NumericVector likelihood = df_nodes["likelihood"];
  vector<double> nodes(likelihood.begin(), likelihood.end());
  vector<vector<unsigned>> edges = make_edges(from, to, names.size());
  Graph g = Graph(nodes, edges);
  vector<unsigned> module = g.random_subgraph(args["module_size"]);
  g.initialize_module(module);
  vector<unsigned> ret = g.onelong_iteration_frequency(args["start"], args["end"]);
  IntegerVector ret_(ret.begin(), ret.end());
  return ret_;
}

bool is_connected(vector<vector<unsigned>> edges, bool from_inner[]) {
  int n = edges.size();
  vector<bool> used(n, false);
  queue<unsigned> q;
  for(int i = 0; i < n; ++i){
    if(from_inner[i]){
      used[i] = true;
      q.push(i);
      break;
    }else if(i == n - 1){
      return true;
    }
  }
  while (!q.empty()) {
    unsigned v = q.front();
    q.pop();
    for (unsigned to : edges[v]) {
      if (from_inner[to] && !used[to]) {
        used[to] = true;
        q.push(to);
      }
    }
  }
  for(int i = 0; i < n; ++i){
    if(from_inner[i] && !used[i]){
      return false;
    }
  }
  return true;
}

// [[Rcpp::export]]
NumericVector real_prob_internal(DataFrame df_edges, DataFrame df_nodes){
  IntegerVector from = df_edges["from"];
  IntegerVector to   = df_edges["to"];
  NumericVector likelihood = df_nodes["likelihood"];
  vector<double> nodes(likelihood.begin(), likelihood.end());
  vector<vector<unsigned>> edges = make_edges(from, to, likelihood.size());
  vector<double> scores(nodes.size(), 0.0);
  double sumscores = 0;
  int n = nodes.size();
  bool x[n];
  for(int i = 1; i < n + 1; ++i){
    for(int j = 0; j < n; ++j)  x[j] = i + j < n ? false : true;
    do{
      if(is_connected(edges, x)){
        double score = 1;
        for(int i = 0; i < n; ++i) if(x[i]) score *= nodes[i];
        sumscores += score;
        for(int i = 0; i < n; ++i) if(x[i]) scores[i] += score;
      }
    } while (next_permutation(x, x + n));
  }
  for(int i = 0; i < n; ++i) scores[i] /= sumscores;
  NumericVector ret(scores.begin(), scores.end());
  return ret;
}
