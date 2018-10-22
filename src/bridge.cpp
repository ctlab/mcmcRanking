#include <Rcpp.h>
#include <vector>
#include <RProgress.h>
#include "mcmc.h"
#include "utils.h"

using namespace Rcpp;
using namespace std;
using mcmc::Graph;

// [[Rcpp::export]]
LogicalVector sample_subgraph_internal(IntegerMatrix edgelist, int gorder, int module_size, size_t niter) {
    RProgress::RProgress pb;
    if (niter > 0) {
        pb = RProgress::RProgress("[:bar] ETA: :eta", niter);
        pb.tick(0);
    }
    vector<double> nodes(gorder, 1);
    Graph g = Graph(nodes, adj_list(edgelist, gorder), true);
    g.initialize_module(g.random_subgraph(module_size));
    for (size_t i = 0; i < niter; ++i) {
        g.next_iteration();
        if (i % 10000 == 9999) {
            Rcpp::checkUserInterrupt();
            pb.tick(10000);
        }
    }
    if (niter > 0)
        pb.tick(niter % 10000);
    LogicalVector ret(gorder, false);
    for (size_t x : g.get_inner_nodes()) {
        ret[x] = true;
    }
    return ret;
}

// [[Rcpp::export]]
NumericVector sample_llh_internal(IntegerMatrix edgelist, NumericVector likelihood, size_t niter, bool fixed_size,
                                  LogicalMatrix start_module) {
    RProgress::RProgress pb;
    if (niter > 0) {
        pb = RProgress::RProgress("[:bar] ETA: :eta", niter);
        pb.tick(0);
    }
    Graph g = Graph(likelihood, adj_list(edgelist, likelihood.size()), fixed_size);
    vector<unsigned> module;
    for (int j = 0; j < start_module.ncol(); ++j) {
        if (start_module(0, j)) {
            module.push_back(j);
        }
    }
    g.initialize_module(module);
    NumericVector llhs(niter, 0);
    for (size_t i = 0; i < niter; ++i) {
        g.next_iteration();
        if (i % 10000 == 9999) {
            Rcpp::checkUserInterrupt();
            pb.tick(10000);
        }
        for (size_t x : g.get_inner_nodes()) {
            llhs[i] += log(likelihood[x]);
        }
    }
    if (niter > 0)
        pb.tick(niter % 10000);
    return llhs;
}

// [[Rcpp::export]]
LogicalMatrix mcmc_sample_internal(IntegerMatrix edgelist, NumericMatrix likelihood, bool fixed_size, size_t niter,
                                   LogicalMatrix start_module) {
    Graph g = Graph((NumericVector) likelihood(_, 0), adj_list(edgelist, likelihood.nrow()), fixed_size);
    size_t order = likelihood.nrow();
    unsigned times = start_module.nrow();
    //LogicalVector ret(order * times, false);
    LogicalMatrix ret(times, order);
    colnames(ret) = rownames(likelihood);

    RProgress::RProgress pb;
    if (times * likelihood.ncol() * niter > 0) {
        pb = RProgress::RProgress("[:bar] ETA: :eta", times * likelihood.ncol() * niter);
        pb.tick(0);
    }
    for (int i = 0; i < times; ++i) {
        vector<unsigned> module;
        for (int j = 0; j < order; ++j) {
            if (start_module(i, j)) {
                module.push_back(j);
            }
        }
        g.initialize_module(module);
        for (int k = 0; k < likelihood.ncol(); ++k) {
            g.set_nodes((NumericVector) likelihood(_, k));
            for (size_t j = 0; j < niter; ++j) {
                g.next_iteration();
                if (j % 10000 == 9999) {
                    Rcpp::checkUserInterrupt();
                    pb.tick(10000);
                }
            }
            if (niter > 0)
                pb.tick(niter % 10000);
        }
        for (size_t x : g.get_inner_nodes()) {
            //ret[x + i * order] = true;
            ret(i, x) = true;
        }
    }
    return ret;
}

// [[Rcpp::export]]
LogicalMatrix
mcmc_onelong_internal(IntegerMatrix edgelist, NumericVector likelihood, bool fixed_size, int module_size, size_t start,
                      size_t niter) {
    RProgress::RProgress pb;
    if (niter > 0) {
        pb = RProgress::RProgress("[:bar] ETA: :eta", niter);
        pb.tick(0);
    }
    Graph g = Graph(likelihood, adj_list(edgelist, likelihood.size()), fixed_size);
    size_t order = likelihood.size();
    g.initialize_module(g.random_subgraph(module_size));
    //LogicalVector ret(order * (niter - start), false);
    LogicalMatrix ret(niter - start, order);
    colnames(ret) = CharacterVector(likelihood.names());
    for (size_t i = 0; i < niter; ++i) {
        g.next_iteration();
        if (i % 10000 == 9999) {
            Rcpp::checkUserInterrupt();
            pb.tick(10000);
        }
        if (i < start) {
            continue;
        }
        for (size_t x : g.get_inner_nodes()) {
            //ret[x + (i - start) * order] = true;
            ret(i-start, x) = true;
        }
    }
    if (niter > 0)
        pb.tick(niter % 10000);
    return ret;
}

// [[Rcpp::export]]
IntegerVector
mcmc_onelong_frequency_internal(IntegerMatrix edgelist, NumericVector likelihood, bool fixed_size, int module_size,
                                size_t start, size_t niter) {
    RProgress::RProgress pb;
    if (niter > 0) {
        pb = RProgress::RProgress("[:bar] ETA: :eta", niter);
        pb.tick(0);
    }
    Graph g = Graph(likelihood, adj_list(edgelist, likelihood.size()), fixed_size);
    size_t order = likelihood.size();
    g.initialize_module(g.random_subgraph(module_size));
    IntegerVector ret(order, 0);
    for (size_t i = 0; i < niter; ++i) {
        g.next_iteration();
        if (i % 10000 == 9999) {
            Rcpp::checkUserInterrupt();
            pb.tick(10000);
        }
        if (i < start) {
            continue;
        }
        for (size_t x : g.get_inner_nodes()) {
            ret[x]++;
        }
    }
    if (niter > 0)
        pb.tick(niter % 10000);
    return ret;
}
