#ifndef MCMC_RANKING_MCMC_H
#define MCMC_RANKING_MCMC_H

#include <random>
#include <vector>

namespace mcmc {
    using namespace std;

    class Graph {
        unsigned number_of_nodes;
        vector<double> nodes;
        vector<vector<unsigned>> edges;

        vector<unsigned> inner_nodes;
        vector<unsigned> outer_nodes;

        bool* from_inner;
        bool* from_outer;

        mt19937 gen;
        uniform_real_distribution<> unirealdis;

        bool update_outer_nodes(unsigned cand_in, unsigned cand_out);

        bool is_connected(unsigned cand_in, unsigned cand_out);

    public:
        Graph(vector<double> nodes, vector<vector<unsigned>> edges);

        //~Graph();

        bool is_connected();

        bool next_iteration(unsigned cand_in, unsigned cand_out, double likelihood_cand_in, double likelihood_cand_out);

        void initialize_module(vector<unsigned> nodes);

        void check_inner_outer();

        vector<unsigned> random_subgraph(unsigned size);

        vector<unsigned> get_inner_nodes();

        vector<unsigned> get_outer_nodes();

        vector<unsigned> subgraph_iteration(unsigned end);

        vector<unsigned> sample_iteration(vector<unsigned> module, unsigned times, unsigned end);

        vector<unsigned> onelong_iteration(unsigned start, unsigned end);

        double onelong_inverse_likelihood(unsigned start, unsigned end);
    };

}

#endif //MCMC_RANKING_MCMC_H
