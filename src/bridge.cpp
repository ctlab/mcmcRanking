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
LogicalVector mcmc_subgraph_internal(DataFrame df_edges, IntegerVector args) {
  IntegerVector from = df_edges["from"];
  IntegerVector to   = df_edges["to"];
  vector<double> nodes(args["nodes_size"], 1);
  vector<vector<unsigned>> edges = make_edges(from, to, args["nodes_size"]);
  Graph g = Graph(nodes, edges, true);
  vector<vector<unsigned>> module;
  module.push_back(g.random_subgraph(args["module_size"]));
  vector<char> ret = g.sample_iteration(module, args["module_size"], 1, args["iter"]);
  LogicalVector ret_(ret.begin(), ret.end());
  return ret_;
}

// [[Rcpp::export]]
NumericVector sample_llh_internal(DataFrame df_edges, DataFrame df_nodes, IntegerVector args, LogicalMatrix start_module) {
  NumericVector likelihood = df_nodes["likelihood"];
  vector<double> nodes(likelihood.begin(), likelihood.end());
  vector<vector<unsigned>> edges = make_edges(df_edges["from"], df_edges["to"], likelihood.size());
  Graph g = Graph(nodes, edges, args["fixed_size"] == 1);
  vector<unsigned> module;
  for(int j = 0; j < start_module.ncol(); ++j){
    if(start_module(0, j)){
      module.push_back(j);
    }
  }
  vector<double> ret = g.sample_llh(module, args["iter"]);
  NumericVector ret_(ret.begin(), ret.end());
  return ret_;
}

// [[Rcpp::export]]
LogicalVector mcmc_sample_internal(DataFrame df_edges, DataFrame df_nodes, IntegerVector args, LogicalMatrix start_module) {
  IntegerVector from = df_edges["from"];
  IntegerVector to   = df_edges["to"];
  IntegerVector names = df_nodes["name"];
  NumericVector likelihood = df_nodes["likelihood"];
  vector<double> nodes(likelihood.begin(), likelihood.end());
  vector<vector<unsigned>> edges = make_edges(from, to, names.size());
  Graph g = Graph(nodes, edges, args["fixed_size"] == 1);
  vector<vector<unsigned>> module;
  for(int i = 0; i < start_module.nrow(); ++i){
    module.push_back(vector<unsigned>());
    for(int j = 0; j < start_module.ncol(); ++j){
      if(start_module(i, j)){
        module[i].push_back(j);
      }
    }
  }
  vector<char> ret = g.sample_iteration(module, args["module_size"], args["times"], args["iter"]);
  LogicalVector ret_(ret.begin(), ret.end());
  return ret_;
}

void dfs(unsigned i, vector<vector<unsigned>> &edges, vector<unsigned> &new_comp, vector<bool> &used){
  used[i] = true;
  new_comp.push_back(i);
  for(unsigned x : edges[i])
    if(!used[x])
      dfs(x, edges, new_comp, used);
}


void cut_points_dfs (int v, int p, int &timer, vector<vector<unsigned>> &edges, vector<unsigned> &ranked, bool used[], int tin[], int fup[], bool is_cut_point[]) {
  used[v] = true;
  tin[v] = timer++;
  fup[v] = timer;
  int children = 0;
  for (size_t i=0; i<edges[v].size(); ++i) {
    int to = edges[v][i];
    if (to == p || ranked[to])
      continue;
    if (used[to]){
      fup[v] = min (fup[v], tin[to]);
    }else {
      cut_points_dfs (to, v, timer, edges, ranked, used, tin, fup, is_cut_point);
      fup[v] = min (fup[v], fup[to]);
      if (fup[to] >= tin[v] && p != -1)
        is_cut_point[v] = true;
      ++children;
    }
  }
  if (p == -1 && children > 1)
    is_cut_point[v] = true;
}


bool get_best_comp(vector<vector<unsigned>> &edges, vector<double> &nodes, vector<unsigned> &ranked, vector<unsigned> &comp, vector<unsigned> &to_remove, int i){
  vector<bool> used(ranked.begin(), ranked.end());
  used[i] = true;
  to_remove.push_back(i);
  for(int j = 0; j < nodes.size(); ++j){
    if(!used[j]){
      vector<unsigned> new_comp;
      dfs(j, edges, new_comp, used);
      if(new_comp.size() > comp.size())
        comp.swap(new_comp);
      for(unsigned x : new_comp){
        to_remove.push_back(x);
        if(nodes[x] > nodes[i])
          return true;
      }
    }
  }
  return false;
}

// [[Rcpp::export]]
IntegerVector probabilistic_rank_internal(DataFrame df_edges, DataFrame df_nodes) {
  NumericVector q = df_nodes["q"];
  unsigned n = q.size();
  vector<double> nodes(q.begin(), q.end());
  vector<vector<unsigned>> edges = make_edges(df_edges["from"], df_edges["to"], n);

  vector<unsigned> ranked(n, 0);
  unsigned not_ranked_n = n;

  double sum_p = 0;
  for(double x : nodes)
    sum_p += 1-x;
  double sum_q = n - sum_p;

  unsigned cand_remove = 0;
  while(not_ranked_n != 0){
    double best_q = sum_q;
    double best_p = -1;

    bool used[n];
    bool is_cut_point[n];
    memset(used, 0, sizeof used);
    memset(is_cut_point, 0, sizeof is_cut_point);
    int tin[n];
    int fup[n];
    int timer = 0;
    for(int i : edges[cand_remove]){
      if(!ranked[i]){
        cut_points_dfs(i, -1, timer, edges, ranked, used, tin, fup, is_cut_point);
        break;
        }
    }
    for(int i = 0; i < n; ++i){
      if(ranked[i])
        continue;
      double cur_p;
      double cur_q;
      if(is_cut_point[i]){
        vector<unsigned> comp;
        vector<unsigned> to_remove;
        if(get_best_comp(edges, nodes, ranked, comp, to_remove, i))
          continue;
        cur_p = 0;
        for(unsigned x : comp)
          cur_p += 1 - nodes[x];
        cur_q = comp.size() - cur_p;
      }else{
        cur_p = sum_p - (1 - nodes[i]);
        cur_q = sum_q - nodes[i];
      }
      double sin = (best_p - sum_p) * (cur_q - sum_q) - (best_q - sum_q) * (cur_p - sum_p);
      bool smaller = cur_p * cur_p + cur_q * cur_q < best_p * best_p + best_q * best_q;
      if(sin > 0 || sin == 0 && smaller || best_p == -1){
        best_p = cur_p;
        best_q = cur_q;
        cand_remove = i;
      }
    }
    if(is_cut_point[cand_remove]){
      vector<unsigned> comp;
      vector<unsigned> to_remove;
      get_best_comp(edges, nodes, ranked, comp, to_remove, cand_remove);
      for(unsigned x : to_remove){
        ranked[x] = not_ranked_n;
      }
      not_ranked_n -= to_remove.size();
    }else{
      ranked[cand_remove] = not_ranked_n--;
    }
    sum_p = best_p;
    sum_q = best_q;
  }
  IntegerVector ret_(ranked.begin(), ranked.end());
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
  Graph g = Graph(nodes, edges, true);
  vector<unsigned> module = g.random_subgraph(args["module_size"]);
  g.initialize_module(module);
  vector<char> ret = g.onelong_iteration(args["start"], args["end"]);
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
  Graph g = Graph(nodes, edges, true);
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
