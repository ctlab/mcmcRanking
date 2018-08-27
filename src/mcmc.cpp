#include <queue>
#include <cmath>
#include "mcmc.h"

namespace mcmc {
    Graph::Graph(vector<double> nodes, vector <vector<unsigned>> edges, bool fixed_size)
            : fixed_size(fixed_size), order(nodes.size()), nodes(nodes), edges(edges), inner(order), outer(order),
              in_nei_c(order, 0), neis(order, unordered_set<size_t>()) {
        random_device rd;
        gen = mt19937(rd());
        unirealdis = uniform_real_distribution<>(0, 1);
    }

    Graph::Graph(Rcpp::NumericVector nodes, vector <vector<unsigned>> edges, bool fixed_size)
            : fixed_size(fixed_size), order(nodes.size()), edges(edges), inner(order), outer(order),
              in_nei_c(order, 0), neis(order, unordered_set<size_t>()) {
        this->nodes = vector<double>(nodes.begin(), nodes.end());
        random_device rd;
        gen = mt19937(rd());
        unirealdis = uniform_real_distribution<>(0, 1);
    }

    void Graph::set_nodes(Rcpp::NumericVector nodes) {
        this->nodes = vector<double>(nodes.begin(), nodes.end());
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
        for (auto &x : neis) {
            x.clear();
        }
        for (unsigned node : nodes) {
            inner.insert(node);
        }
        for (size_t i = 0; i < inner.size(); ++i) {
            unsigned v = inner.get(i);
            for (unsigned neighbour : edges[v]) {
                in_nei_c[neighbour]++;
                if (inner.contains(neighbour)) {
                    neis[v].insert(neighbour);
                } else if (!outer.contains(neighbour)) {
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
            for (auto &to : neis[v]) {
                unsigned u = inner.get_index(to);
                if (u != -1) {
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

    void Graph::inner_update(unsigned v, bool is_erased) {
        if (is_erased) {
            inner.erase(v);
            for (auto &neighbour : neis[v]) {
                neis[neighbour].erase(v);
            }
            neis[v].clear();
        } else {
            inner.insert(v);
            for (unsigned neighbour : edges[v]) {
                if (inner.contains(neighbour)) {
                    neis[v].insert(neighbour);
                    neis[neighbour].insert(v);
                }
            }
        }
    }

    void Graph::update_neighbours(unsigned v, bool is_erased) {
        if (is_erased) {
            if (inner.size() != 0)
                outer.insert(v);
            for (unsigned neighbour : edges[v]) {
                if (--in_nei_c[neighbour] == 0) {
                    if (!inner.contains(neighbour)) {
                        outer.erase(neighbour);
                    }
                }
            }
        } else {
            if (inner.size() != 1)
                outer.erase(v);
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
                       (outer.size() <= edges[cand_in].size() - in_nei_c[cand_in] ?
                        1 : outer.size() - edges[cand_in].size() + in_nei_c[cand_in]);
            if (gen_p >= p)
                return false;
            inner_update(cand_in, true);
            inner_update(cand_out, false);
            if (!is_connected()) {
                inner_update(cand_out, true);
                inner_update(cand_in, false);
                return false;
            }
            unsigned cur_size_outer = outer.size();
            update_outer_nodes(cand_in, cand_out);
            unsigned new_size_outer = outer.size();
            p = (nodes[cand_out] * cur_size_outer) / (nodes[cand_in] * new_size_outer);
            if (gen_p < p) {
                return true;
            }
            inner_update(cand_out, true);
            inner_update(cand_in, false);
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
                int new_in_out_size = inner.size() + outer.size() - edges[cand].size() + in_nei_c[cand];
                if (inner.size() == 1) {
                    new_in_out_size = order;
                } else if (cur_in_out_size <= edges[cand].size() - in_nei_c[cand]) {
                    new_in_out_size = 1;
                }
                double p = cur_in_out_size / (nodes[cand] * new_in_out_size);
                if (gen_p >= p) {
                    return false;
                }
                inner_update(cand, erase);
                if (!is_connected()) {
                    inner_update(cand, !erase);
                    return false;
                }
            } else {
                double p = nodes[cand] * (inner.size() == 0 ? 1.0 * order / (1 + edges[cand].size()) : 1);
                if (gen_p >= p) {
                    return false;
                }
                inner_update(cand, erase);
            }
            update_neighbours(cand, erase);
            unsigned new_in_out_size = inner.size() == 0 ? order : outer.size() + inner.size();
            double p = (erase ? 1 / nodes[cand] : nodes[cand]) * (1.0 * cur_in_out_size / new_in_out_size);
            if (gen_p < p) {
                return true;
            }
            inner_update(cand, !erase);
            update_neighbours(cand, !erase);
            return false;
        }
    }

    vector <size_t> Graph::get_inner_nodes() {
        return inner.get_all();
    }

    vector <size_t> Graph::get_outer_nodes() {
        return outer.get_all();
    }
};
