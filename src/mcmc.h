#ifndef MCMC_RANKING_MCMC_H
#define MCMC_RANKING_MCMC_H

#include <random>
#include <vector>
#include <unordered_map>
#include "hsa.h"

namespace mcmc {
    using namespace std;
    class Graph {
        size_t order;
        vector<double> nodes;
        vector<vector<unsigned>> edges;

        HSA inner;
        HSA outer;

        mt19937 gen;
        uniform_real_distribution<> unirealdis;

        bool update_outer_nodes(unsigned cand_in, unsigned cand_out);

    public:
        Graph(vector<double> nodes, vector<vector<unsigned>> edges);

        bool is_connected();

        bool next_iteration();

        void initialize_module(vector<unsigned> nodes);

        vector<unsigned> random_subgraph(size_t size);

        vector<unsigned> get_inner_nodes();

        vector<unsigned> get_outer_nodes();

        vector<unsigned> sample_iteration(vector<unsigned> module, size_t times, size_t end);

        vector<unsigned> onelong_iteration(size_t start, size_t end);

        vector<unsigned> onelong_iteration_frequency(size_t start, size_t end);
    };

}

#endif //MCMC_RANKING_MCMC_H
