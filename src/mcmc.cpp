#include <queue>
#include <unordered_set>
#include <cmath>
#include "mcmc.h"
#include <Rcpp.h>

namespace mcmc {
    Graph::Graph(vector<double> nodes, vector<vector<unsigned>> edges, bool fixed_size)
            : fixed_size(fixed_size), order(nodes.size()), nodes(nodes), edges(edges), inner(order), outer(order),
              in_nei_c(order, 0) {
        random_device rd;
        gen = mt19937(rd());
        unirealdis = uniform_real_distribution<>(0, 1);
    }

    Graph::Graph(Rcpp::NumericVector nodes, vector<vector<unsigned>> edges, bool fixed_size)
            : fixed_size(fixed_size), order(nodes.size()), edges(edges), inner(order), outer(order),
              in_nei_c(order, 0) {
        this->nodes = vector<double>(nodes.begin(), nodes.end());
        random_device rd;
        gen = mt19937(rd());
        unirealdis = uniform_real_distribution<>(0, 1);
    }

    vector<unsigned> Graph::random_subgraph(size_t size) {
        if (size == 0) {
            return vector<unsigned>();
        }
        unordered_set<unsigned> sg;
        HSA candidates(order);
        size_t ind = uniform_int_distribution<>(0, order - 1)(gen);
        sg.insert(ind);
        for (unsigned cand_node : edges[ind]) {
            if (sg.count(cand_node) < 1)
                candidates.insert(cand_node);
        }
        while (sg.size() != size) {
            ind = uniform_int_distribution<>(0, candidates.size() - 1)(gen);
            unsigned new_node = candidates.get(ind);
            for (unsigned neigbour : edges[new_node]) {
                if (!candidates.contains(neigbour)) {
                    if (sg.count(neigbour) < 1) {
                        candidates.insert(neigbour);
                    }
                }
            }
            sg.insert(new_node);
            candidates.erase(new_node);
        }
        return vector<unsigned>(sg.begin(), sg.end());
    }

    void Graph::initialize_module(vector<unsigned> nodes) {
        inner.clear();
        outer.clear();
        std::fill(in_nei_c.begin(), in_nei_c.end(), 0);
        for (unsigned node : nodes) {
            inner.insert(node);
        }
        for (size_t i = 0; i < inner.size(); ++i) {
            for (unsigned neighbour : edges[inner.get(i)]) {
                in_nei_c[neighbour]++;
                if ((!inner.contains(neighbour)) && (!outer.contains(neighbour))) {
                    outer.insert(neighbour);
                }
            }
        }
    }

    bool Graph::is_connected() {
        if (inner.size() == 0)
            return true;
        vector<char> used(inner.size(), false);
        queue<unsigned> q;
        int el = inner.get(0);
        used[0] = true;
        q.push(el);
        while (!q.empty()) {
            unsigned v = q.front();
            q.pop();
            for (unsigned to : edges[v]) {
                unsigned u = inner.get_index(to);
                if(u != -1) {
                    if (!used[u]) {
                        used[u] = true;
                        q.push(to);
                    }
                }
            }
        }
        for (size_t i = 0; i < inner.size(); ++i) {
            if (!used[i]) {
                return false;
            }
        }
        return true;
    }

    void Graph::update_outer_nodes(unsigned cand_in, unsigned cand_out) {
        outer.swap(cand_out, cand_in);

        for (unsigned neighbour : edges[cand_out]) {
            if (in_nei_c[neighbour]++ == 0 && neighbour != cand_in) {
                if (!inner.contains(neighbour)) {
                    outer.insert(neighbour);
                }
            }
        }
        for (unsigned neighbour : edges[cand_in]) {
            if (--in_nei_c[neighbour] == 0 && neighbour != cand_out) {
                if (!inner.contains(neighbour)) {
                    outer.erase(neighbour);
                }
            }
        }
    }

    void Graph::update_neighbours(unsigned v, bool is_erased) {
        if (is_erased) {
            for (unsigned neighbour : edges[v]) {
                if (--in_nei_c[neighbour] == 0) {
                    if (!inner.contains(neighbour)) {
                        outer.erase(neighbour);
                    }
                }
            }
        } else {
            for (unsigned neighbour : edges[v]) {
                if (in_nei_c[neighbour]++ == 0) {
                    if (!inner.contains(neighbour)) {
                        outer.insert(neighbour);
                    }
                }
            }
        }
    }

    bool Graph::next_iteration() {
        if (fixed_size) {
            if (inner.size() == 0) {
                return true;
            }
            if (inner.size() == order) {
                return false;
            }
            unsigned cand_in = inner.get(uniform_int_distribution<>(0, inner.size() - 1)(gen));
            unsigned cand_out = outer.get(uniform_int_distribution<>(0, outer.size() - 1)(gen));
            double gen_p = unirealdis(gen);
            double p = (nodes[cand_out] * outer.size()) / nodes[cand_in] *
                       (outer.size() < edges[cand_in].size() ? 1 : outer.size() - edges[cand_in].size() + 1);
            if (gen_p >= p)
                return false;
            inner.swap(cand_in, cand_out);
            if (!is_connected()) {
                inner.swap(cand_out, cand_in);
                return false;
            }
            unsigned cur_size_outer = outer.size();
            update_outer_nodes(cand_in, cand_out);
            unsigned new_size_outer = outer.size();
            p = (nodes[cand_out] * cur_size_outer) / (nodes[cand_in] * new_size_outer);
            if (gen_p < p) {
                return true;
            }
            inner.swap(cand_out, cand_in);
            update_outer_nodes(cand_out, cand_in);
            return false;
        } else {
            unsigned cur_in_out_size = inner.size() == 0 ? order : outer.size() + inner.size();
            bool erase = unirealdis(gen) < (1.0 * inner.size()) / cur_in_out_size;
            unsigned cand;
            if (erase) {
                cand = inner.get(uniform_int_distribution<>(0, inner.size() - 1)(gen));
            } else {
                cand = outer.size() == 0 ?
                       uniform_int_distribution<>(0, order - 1)(gen) :
                       outer.get(uniform_int_distribution<>(0, outer.size() - 1)(gen));
            }
            double gen_p = unirealdis(gen);
            if (erase) {
                int new_in_out_size = inner.size() + outer.size() - edges[cand].size();
                if (inner.size() == 0) {
                    new_in_out_size = order;
                } else if (inner.size() + outer.size() <= edges[cand].size()) {
                    new_in_out_size = 1;
                }
                double p = (inner.size() + outer.size()) / (nodes[cand] * new_in_out_size);
                if (gen_p >= p) {
                    return false;
                }
                inner.erase(cand);
                if (!is_connected()) {
                    inner.insert(cand);
                    return false;
                }
                if (inner.size() != 0)
                    outer.insert(cand);
            } else {
                double p = nodes[cand] * (inner.size() == 0 ? 1.0 * order / (1 + edges[cand].size()) : 1);
                if (gen_p >= p) {
                    return false;
                }
                if (inner.size() != 0)
                    outer.erase(cand);
                inner.insert(cand);
            }
            update_neighbours(cand, erase);
            unsigned new_in_out_size = inner.size() == 0 ? order : outer.size() + inner.size();
            double p = (erase ? 1 / nodes[cand] : nodes[cand]) * (1.0 * cur_in_out_size / new_in_out_size);
            if (gen_p < p) {
                return true;
            }
            if (erase) {
                inner.insert(cand);
                if (inner.size() != 1)
                    outer.erase(cand);
            } else {
                if (inner.size() != 1)
                    outer.insert(cand);
                inner.erase(cand);
            }
            update_neighbours(cand, !erase);
            return false;
        }
    }

    vector<size_t> Graph::get_inner_nodes() {
        return inner.get_all();
    }

    vector<size_t> Graph::get_outer_nodes() {
        return outer.get_all();
    }

    vector<double> Graph::sample_llh(vector<unsigned> module, size_t end) {
        vector<double> likelihoods(end);
        initialize_module(module);
        for (size_t i = 0; i < end; ++i) {
            next_iteration();
            likelihoods[i] = 0;
            for (unsigned x : inner.get_all()) {
                likelihoods[i] += log(nodes[x]);
            }
        }
        return likelihoods;
    }

    vector<char> Graph::sample_iteration(vector<vector<unsigned>> module, size_t times, size_t end) {
        vector<char> ret(order * times, false);
        for (size_t i = 0; i < times; ++i) {
            Rcpp::checkUserInterrupt();
            initialize_module(module[i]);
            for (size_t j = 0; j < end; ++j) {
                next_iteration();
            }
            for (unsigned x : inner.get_all()) {
                ret[x + i * order] = true;
            }
        }
        return ret;
    }

    vector<char> Graph::onelong_iteration(size_t start, size_t end) {
        vector<char> ret(order * (end - start), false);
        for (size_t i = 0; i < end; ++i) {
            next_iteration();
            if (i >= start) {
                for (unsigned x : inner.get_all()) {
                    ret[x + (i - start) * order] = true;
                }
            }
        }
        return ret;
    }

    vector<unsigned> Graph::onelong_iteration_frequency(size_t start, size_t end) {
        vector<unsigned> ret(order, 0);
        for (size_t i = 0; i < end; ++i) {
            next_iteration();
            if (i >= start) {
                for (unsigned x : inner.get_all()) {
                    ret[x]++;
                }
            }
        }
        return ret;
    }
};
