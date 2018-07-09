#ifndef MCMC_RANKING_MCMC_H
#define MCMC_RANKING_MCMC_H

#include <random>
#include <vector>
#include <unordered_map>
#include "hsa.h"
#include <Rcpp.h>

namespace mcmc {
    using namespace std;
    class Graph {
        bool fixed_size;
        size_t order;
        vector<double> nodes;
        vector<vector<unsigned>> edges;

        HSA inner;
        HSA outer;

        mt19937 gen;
        uniform_real_distribution<> unirealdis;

        void update_outer_nodes(unsigned cand_in, unsigned cand_out);

        void update_neighbours(unsigned v, bool is_erased);

    public:
        Graph(vector<double> nodes, vector<vector<unsigned>> edges, bool fixed_size);

        Graph(Rcpp::NumericVector nodes, vector<vector<unsigned>> edges, bool fixed_size);

        bool is_connected();

        bool next_iteration();

        void initialize_module(vector<unsigned> nodes);

        vector<unsigned> random_subgraph(size_t size);

        vector<unsigned> get_inner_nodes();

        vector<unsigned> get_outer_nodes();

        vector<double> sample_llh(vector<unsigned> module, size_t end);

        vector<char> sample_iteration(vector<vector<unsigned>> module, size_t module_size, size_t end);

        vector<char> onelong_iteration(size_t start, size_t end);

        vector<unsigned> onelong_iteration_frequency(size_t start, size_t end);
    };

}

#endif //MCMC_RANKING_MCMC_H
