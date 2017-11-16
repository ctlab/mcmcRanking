#include <Rcpp.h>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include "mcmc-ranking/src/mcmc.h"
using namespace Rcpp;
using namespace std;
using mcmc::Graph;
// This is a simple function using Rcpp that creates an R list
// containing a character vector and a numeric vector.
//
// Learn more about how to use Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//
// and browse examples of code using Rcpp at:
//
//   http://gallery.rcpp.org/
//

pair<map<string, unsigned>, map<unsigned,string>> make_bimap(CharacterVector names){
  map<string, unsigned> name_to_id;
  for(int i = 0; i < names.size(); ++i){
    name_to_id.insert(pair<string, unsigned>(as<string>(names[i]), i));
  }
  map<unsigned, string> id_to_name;
  for (map<string, unsigned>::iterator i = name_to_id.begin(); i != name_to_id.end(); ++i)
    id_to_name[i->second] = i->first;
  return pair<map<string, unsigned>, map<unsigned,string>>(name_to_id, id_to_name);
}

// [[Rcpp::export]]
CharacterVector mcmc_subgraph(DataFrame df_edges, IntegerVector size, IntegerVector iter) {
  CharacterVector from = df_edges["from"];
  CharacterVector to   = df_edges["to"];
  CharacterVector nodes_ = union_(from, to);
  pair<map<string, unsigned>, map<unsigned,string>> bimap = make_bimap(nodes_);

  vector<double> nodes(nodes_.size());
  vector<vector<unsigned>> edges(nodes_.size());
  for(int i = 0; i < from.size(); ++i){
    unsigned from_ = bimap.first.find(as<string>(from[i])) -> second;
    unsigned to_   = bimap.first.find(as<string>(to[i])) -> second;
    edges[from_].push_back(to_);
    edges[to_].push_back(from_);
  }
  Graph g = Graph(nodes, edges);
  g.initialize_module(g.random_subgraph(size[0]));
  vector<unsigned> ret = g.subgraph_iteration(iter[0]);
  CharacterVector ret_;
  for(int i = 0; i < ret.size(); ++i){
    ret_.push_back(bimap.second.find(ret[i])->second);
  }
  return ret_;
}


// [[Rcpp::export]]
NumericVector mcmc_sample(DataFrame df_edges, DataFrame df_nodes, IntegerVector size, IntegerVector times, IntegerVector iter) {
  CharacterVector from = df_edges["from"];
  CharacterVector to   = df_edges["to"];
  CharacterVector names = df_nodes["name"];
  NumericVector likelihood = df_nodes["likelihood"];
  pair<map<string, unsigned>, map<unsigned,string>> bimap = make_bimap(names);

  vector<double> nodes(names.size());
  for(int i = 0; i < names.size(); ++i){
    nodes[i] = likelihood[i];
  }
  vector<vector<unsigned>> edges(names.size());
  for(int i = 0; i < from.size(); ++i){
    unsigned from_ = bimap.first.find(as<string>(from[i])) -> second;
    unsigned to_   = bimap.first.find(as<string>(to[i])) -> second;
    edges[from_].push_back(to_);
    edges[to_].push_back(from_);
  }
  Graph g = Graph(nodes, edges);
  vector<unsigned> module = g.random_subgraph(size[0]);
  vector<double> ret = g.sample_iteration(module, times[0], iter[0]);
  NumericVector ret_;
  for(int i = 0; i < ret.size(); ++i){
    ret_.push_back(ret[i]);
  }
  ret_.names() = names;
  return ret_;
}


// [[Rcpp::export]]
NumericVector mcmc_onelong(DataFrame df_edges, DataFrame df_nodes, IntegerVector size, IntegerVector start, IntegerVector end) {
  CharacterVector from = df_edges["from"];
  CharacterVector to   = df_edges["to"];
  CharacterVector names = df_nodes["name"];
  NumericVector likelihood = df_nodes["likelihood"];
  pair<map<string, unsigned>, map<unsigned,string>> bimap = make_bimap(names);

  vector<double> nodes(names.size());
  for(int i = 0; i < names.size(); ++i){
    nodes[i] = likelihood[i];
  }
  vector<vector<unsigned>> edges(names.size());
  for(int i = 0; i < from.size(); ++i){
    unsigned from_ = bimap.first.find(as<string>(from[i])) -> second;
    unsigned to_   = bimap.first.find(as<string>(to[i])) -> second;
    edges[from_].push_back(to_);
    edges[to_].push_back(from_);
  }
  Graph g = Graph(nodes, edges);
  vector<unsigned> module = g.random_subgraph(size[0]);
  g.initialize_module(module);
  vector<double> ret = g.onelong_iteration(start[0], end[0]);
  NumericVector ret_;
  for(int i = 0; i < ret.size(); ++i){
    ret_.push_back(ret[i]);
  }
  ret_.names() = names;
  return ret_;
}


// [[Rcpp::export]]
NumericVector mcmc_reciprocal_likelihood(DataFrame df_edges, DataFrame df_nodes, IntegerVector size, IntegerVector start, IntegerVector end) {
  CharacterVector from = df_edges["from"];
  CharacterVector to   = df_edges["to"];
  CharacterVector names = df_nodes["name"];
  NumericVector likelihood = df_nodes["likelihood"];
  pair<map<string, unsigned>, map<unsigned,string>> bimap = make_bimap(names);

  vector<double> nodes(names.size());
  for(int i = 0; i < names.size(); ++i){
    nodes[i] = likelihood[i];
  }
  vector<vector<unsigned>> edges(names.size());
  for(int i = 0; i < from.size(); ++i){
    unsigned from_ = bimap.first.find(as<string>(from[i])) -> second;
    unsigned to_   = bimap.first.find(as<string>(to[i])) -> second;
    edges[from_].push_back(to_);
    edges[to_].push_back(from_);
  }
  Graph g = Graph(nodes, edges);
  vector<unsigned> module = g.random_subgraph(size[0]);
  g.initialize_module(module);
  return(NumericVector::create(g.onelong_reciprocal_likelihood(start[0], end[0])));
}
