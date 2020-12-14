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

        // neis[v][i].first - neighpour
        // neis[v][i].second - index of v in neis[neighbour]
        vector <vector<pair<size_t, size_t>>> neis;

        mt19937 gen;
        uniform_real_distribution<> unirealdis;

        // These vectors are not local in is_connected in order to avoid memory allocations
        vector<pair<unsigned, unsigned>> bfsUsed;
        vector<unsigned> dsu, dsuCnt;
        unsigned bfsUsedIteration = 0;
        vector<size_t> bfsQueue;

        void update_outer_nodes(unsigned cand_in, unsigned cand_out);

        void update_neighbours(unsigned v, bool is_erased);

        void inner_update(unsigned v, bool is_erased);

        bool is_connected(size_t erased);

    public:
        Graph(vector<double> nodes, vector <vector<unsigned>> edges, bool fixed_size);

        Graph(Rcpp::NumericVector nodes, vector <vector<unsigned>> edges, bool fixed_size);

        void set_nodes(Rcpp::NumericVector nodes);

        bool next_iteration();

        void initialize_module(vector<unsigned> nodes);

        vector<unsigned> random_subgraph(size_t size);

        vector <size_t> get_inner_nodes();

        vector <size_t> get_outer_nodes();
    };

}

#endif //MCMC_RANKING_MCMC_H
