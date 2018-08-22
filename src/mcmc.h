#ifndef MCMC_RANKING_MCMC_H
#define MCMC_RANKING_MCMC_H

#include <random>
#include <vector>
#include <unordered_set>
#include <Rcpp.h>
#include "hsa.h"

namespace mcmc {
    using namespace std;

    class Graph {
        bool fixed_size;
        size_t order;
        vector<double> nodes;
        vector <vector<unsigned>> edges;

        HSA inner;
        HSA outer;

        vector <size_t> in_nei_c;
        vector <unordered_set<size_t>> neis;

        mt19937 gen;
        uniform_real_distribution<> unirealdis;

        void update_outer_nodes(unsigned cand_in, unsigned cand_out);

        void update_neighbours(unsigned v, bool is_erased);

        void inner_update(unsigned v, bool is_erased);

    public:
        Graph(vector<double> nodes, vector <vector<unsigned>> edges, bool fixed_size);

        Graph(Rcpp::NumericVector nodes, vector <vector<unsigned>> edges, bool fixed_size);

        bool is_connected();

        bool next_iteration();

        void initialize_module(vector<unsigned> nodes);

        vector<unsigned> random_subgraph(size_t size);

        vector <size_t> get_inner_nodes();

        vector <size_t> get_outer_nodes();
    };

}

#endif //MCMC_RANKING_MCMC_H
