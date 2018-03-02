#include <queue>
#include <iostream>
#include <algorithm>
#include "mcmc.h"

// TODO: deconstructor
namespace mcmc {
    Graph::Graph(vector<double> nodes, vector<vector<unsigned>> edges) {
        this->number_of_nodes = nodes.size();
        this->nodes = nodes;
        this->edges = edges;
        random_device rd;
        gen = mt19937(rd());
        unirealdis = uniform_real_distribution<>(0, 1);
    }

    vector<unsigned> Graph::random_subgraph(unsigned size) {
        vector<unsigned> sg;
        vector<unsigned> candidates;
        unsigned ind = uniform_int_distribution<>(0, number_of_nodes - 1)(gen);
        sg.push_back(ind);
        for (unsigned cand_node : edges[ind]) {
            if(cand_node != ind)
                candidates.push_back(cand_node);
        }
        while (sg.size() != size) {
            ind = candidates.size() == 1 ? 0 : uniform_int_distribution<>(0, candidates.size() - 1)(gen);
            unsigned new_node = candidates[ind];
            for (unsigned neigbour : edges[new_node]) {
                for (int i = 0; i < candidates.size(); ++i) {
                    if (neigbour == candidates[i]) {
                        break;
                    } else if (i == candidates.size() - 1) {
                        for (int j = 0; j < sg.size(); ++j) {
                            if (neigbour == sg[j])
                                break;
                            else if (j == sg.size() - 1)
                                candidates.push_back(neigbour);
                        }
                    }
                }
            }
            sg.push_back(new_node);
            candidates.erase(candidates.begin() + ind);
        }
        return sg;
    }

    void Graph::initialize_module(vector<unsigned> nodes) {
        from_inner = new bool[number_of_nodes]{false};
        from_outer = new bool[number_of_nodes]{false};
        inner_nodes.clear();
        outer_nodes.clear();
        for (unsigned node : nodes) {
            inner_nodes.push_back(node);
            from_inner[node] = true;
        }
        for (unsigned node : inner_nodes) {
            for (unsigned neighbour : edges[node]) {
                if (!from_inner[neighbour] && !from_outer[neighbour]) {
                    outer_nodes.push_back(neighbour);
                    from_outer[neighbour] = true;
                }
            }
        }
    }

    bool Graph::is_connected() {
        vector<bool> used(number_of_nodes, false);
        queue<unsigned> q;
        used[inner_nodes[0]] = true;
        q.push(inner_nodes[0]);
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
        for (unsigned node : inner_nodes) {
            if (!used[node]) {
                return false;
            }
        }
        return true;
    }

    bool Graph::is_connected(unsigned cand_in, unsigned cand_out) {
        from_inner[cand_in] = false;
        from_inner[cand_out] = true;

        vector<bool> used(number_of_nodes, false);
        queue<unsigned> q;
        used[cand_out] = true;
        q.push(cand_out);
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

        from_inner[cand_in] = true;
        from_inner[cand_out] = false;
        for (unsigned node : inner_nodes) {
            if (!used[node] && node != cand_in) {
                return false;
            }
        }
        if (!used[cand_out]) {
            return false;
        }
        return true;
    }

    bool Graph::update_outer_nodes(unsigned cand_in, unsigned cand_out) {
        for (int i = 0; i < inner_nodes.size(); ++i) {
            if (inner_nodes[i] == cand_in) {
                from_inner[cand_in] = false;
                from_inner[cand_out] = true;
                inner_nodes[i] = cand_out;
            }
        }

        for (int i = 0; i < outer_nodes.size(); ++i) {
            if (outer_nodes[i] == cand_out) {
                from_outer[cand_in] = true;
                from_outer[cand_out] = false;
                outer_nodes[i] = cand_in;
            }
        }

        for (unsigned neighbour : edges[cand_out]) {
            if (!from_inner[neighbour] && !from_outer[neighbour]) {
                from_outer[neighbour] = true;
                outer_nodes.push_back(neighbour);
            }
        }

        for (unsigned neighbour : edges[cand_in]) {
            if (from_inner[neighbour]) {
                continue;
            }
            for (int j = 0; j < edges[neighbour].size(); ++j) {
                if (from_inner[edges[neighbour][j]]) {
                    break;
                } else if (j == edges[neighbour].size() - 1) {
                    for (int k = 0; k < outer_nodes.size(); ++k) {
                        if (outer_nodes[k] == neighbour) {
                            outer_nodes.erase(outer_nodes.begin() + k);
                            from_outer[neighbour] = false;
                            break;
                        }
                    }
                }
            }
        }
    }

    bool Graph::next_iteration(unsigned cand_in, unsigned cand_out, double likelihood_cand_in,
                               double likelihood_cand_out) {
        if (!is_connected(cand_in, cand_out)) {
            return false;
        }
        unsigned cur_size_outer = outer_nodes.size();
        update_outer_nodes(cand_in, cand_out);
        unsigned new_size_outer = outer_nodes.size();
        double p = (likelihood_cand_out * cur_size_outer) / (likelihood_cand_in * new_size_outer);
        if (unirealdis(gen) < p) {
            return true;
        }
        update_outer_nodes(cand_out, cand_in);
        return false;
    }


    void Graph::check_inner_outer() {
        for (int i = 0; i < inner_nodes.size(); ++i) {
            for (int j = i + 1; j < inner_nodes.size(); ++j) {
                if (inner_nodes[i] == inner_nodes[j])
                    std::cout << "inner_nodes has doublicate elements!\n";
            }
        }
        for (int i = 0; i < outer_nodes.size(); ++i) {
            for (int j = i + 1; j < outer_nodes.size(); ++j) {
                if (outer_nodes[i] == outer_nodes[j])
                    std::cout << "outer_nodes has doublicate elements!\n";
            }
        }
        for (int i = 0; i < inner_nodes.size(); ++i) {
            for (int j = 0; j < outer_nodes.size(); ++j) {
                if (inner_nodes[i] == outer_nodes[j])
                    std::cout << "Node belongs to inner_nodes and outerNOdes.!\n";
            }
        }
        if (!is_connected())
            std::cout << "The graph is not connected.\n";
    }

    vector<unsigned> Graph::get_inner_nodes() {
        return inner_nodes;
    }

    vector<unsigned> Graph::get_outer_nodes() {
        return outer_nodes;
    }

    vector<unsigned> Graph::subgraph_iteration(unsigned end) {
        for (int i = 0; i < end; ++i) {
            unsigned candidate_in = inner_nodes[uniform_int_distribution<>(0, inner_nodes.size() - 1)(gen)];
            unsigned candidate_out = outer_nodes[uniform_int_distribution<>(0, outer_nodes.size() - 1)(gen)];
            next_iteration(candidate_in, candidate_out, 1.0, 1.0);
        }
        return inner_nodes;
    }

    vector<unsigned> Graph::sample_iteration(vector<unsigned> module, unsigned times, unsigned end) {
        vector<unsigned> ret;
        for (int i = 0; i < times; ++i) {
            initialize_module(module);
            for (int j = 0; j < end; ++j) {
                unsigned candidate_in = inner_nodes[uniform_int_distribution<>(0, inner_nodes.size() - 1)(gen)];
                unsigned candidate_out = outer_nodes[uniform_int_distribution<>(0, outer_nodes.size() - 1)(gen)];
                next_iteration(candidate_in, candidate_out, nodes[candidate_in], nodes[candidate_out]);
            }
            for (unsigned x : inner_nodes) {
                ret.push_back(x);
            }
        }
        return ret;
    }

    vector<unsigned> Graph::onelong_iteration(unsigned start, unsigned end) {
        vector<unsigned> ret;
        for (int i = 0; i < end; ++i) {
            unsigned candidate_in = inner_nodes[uniform_int_distribution<>(0, inner_nodes.size() - 1)(gen)];
            unsigned candidate_out = outer_nodes[uniform_int_distribution<>(0, outer_nodes.size() - 1)(gen)];
            next_iteration(candidate_in, candidate_out, nodes[candidate_in], nodes[candidate_out]);
            if (i >= start) {
                for (unsigned x : inner_nodes) {
                    ret.push_back(x);
                }
            }
        }
        return ret;
    }

    double Graph::onelong_inverse_likelihood(unsigned start, unsigned end) {
        vector<double> likelihood;
        double cur_likelihood = 1;
        for (unsigned i : inner_nodes) {
            cur_likelihood *= nodes[i];
        }
        for (int i = 0; i < end; ++i) {
            unsigned candidate_in = inner_nodes[uniform_int_distribution<>(0, inner_nodes.size() - 1)(gen)];
            unsigned candidate_out = outer_nodes[uniform_int_distribution<>(0, outer_nodes.size() - 1)(gen)];
            if (next_iteration(candidate_in, candidate_out, nodes[candidate_in], nodes[candidate_out])) {
                cur_likelihood *= nodes[candidate_out] / nodes[candidate_in];
                if (i >= start) {
                    likelihood.push_back(cur_likelihood);
                }
            }
        }
        std::sort(likelihood.begin(), likelihood.end());
        double mean = 0;
        for (int i = likelihood.size() - 1; i >= 0; --i) {
            mean += 1 / likelihood[i];
        }
        mean /= end - start;
        return mean;
    }

};
