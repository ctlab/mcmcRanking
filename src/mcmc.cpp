#include <queue>
#include <cmath>
#include "mcmc.h"

namespace mcmc {
    Graph::Graph(vector<double> nodes, vector <vector<unsigned>> edges, bool fixed_size)
            : fixed_size(fixed_size), order(nodes.size()), nodes(nodes), edges(edges), inner(order), outer(order),
              in_nei_c(order, 0), neis(order),
              bfsUsed(nodes.size(), {0, 0}) {
        random_device rd;
        gen = mt19937(rd());
        unirealdis = uniform_real_distribution<>(0, 1);
        bfsQueue.reserve(nodes.size());
    }

    Graph::Graph(Rcpp::NumericVector nodes, vector <vector<unsigned>> edges, bool fixed_size)
            : Graph(vector<double>(nodes.begin(), nodes.end()), edges, fixed_size) {}

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
                    if (v < neighbour) {
                        neis[v].emplace_back(neighbour, neis[neighbour].size());
                        neis[neighbour].emplace_back(v, neis[v].size() - 1);
                    }
                } else if (!outer.contains(neighbour)) {
                    outer.insert(neighbour);
                }
            }
        }
    }

    static unsigned dsuGetRoot(vector<unsigned> & dsu, unsigned v) {
        auto &p = dsu[v];
        return p == v ? v : p = dsuGetRoot(dsu, p);
    }

    bool Graph::is_connected(size_t erased) {
        size_t k = neis[erased].size();
        if (k <= 1) return true;
        for (auto const& x : neis[erased]) if (neis[x.first].size() == 1) return false;
        ++bfsUsedIteration;
        bfsQueue.clear();
        for (size_t i = 0; i < k; ++i) {
            size_t v = neis[erased][i].first;
            bfsUsed[v] = {bfsUsedIteration, i};
            bfsQueue.push_back(v);
        }
        dsu.clear();
        dsu.resize(k);
        iota(dsu.begin(), dsu.end(), 0);
        dsuCnt.clear();
        dsuCnt.resize(k, 1);

        size_t compCount = dsu.size();
        for (size_t i = 0; i < bfsQueue.size(); ++i) {
            auto v = bfsQueue[i];
            auto comp = dsuGetRoot(dsu, bfsUsed[v].second);
            auto curDsuCnt = dsuCnt[comp] - 1;
            for (auto const& tov : neis[v]) {
                auto to = tov.first;
                if (to == erased) continue;
                auto &used = bfsUsed[to];
                if (used.first != bfsUsedIteration) {
                    if (neis[to].size() > 1) {
                        used = {bfsUsedIteration, comp};
                        bfsQueue.push_back(to);
                        ++curDsuCnt;
                    }
                } else {
                    auto comp1 = dsuGetRoot(dsu, used.second);
                    if (comp1 != comp) {
                        --compCount;
                        if (compCount == 1) return true;
                        dsu[comp1] = comp;
                        curDsuCnt += dsuCnt[comp1];
                    }
                }
            }
            if (curDsuCnt == 0) return false;
            dsuCnt[comp] = curDsuCnt;
        }
        return false;
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
            for (auto const& v2p : neis[v]) {
                auto v2 = v2p.first;
                auto pos = v2p.second;
                auto &vec = neis[v2];
                auto newP = vec[pos] = neis[v2].back();
                neis[newP.first][newP.second].second = pos;
                vec.pop_back();
            }
            neis[v].clear();
        } else {
            inner.insert(v);
            for (unsigned v2 : edges[v]) {
                if (inner.contains(v2)) {
                    neis[v].emplace_back(v2, neis[v2].size());
                    neis[v2].emplace_back(v, neis[v].size() - 1);
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
            unsigned cur_size_outer = outer.size();
            unsigned new_size_outer;
            if (inner.size() == 1) {
                new_size_outer = edges[cand_out].size();
            } else {
                new_size_outer = outer.size();
                for (auto x : edges[cand_in]) {
                    if (!--in_nei_c[x]) --new_size_outer;
                }
                for (auto x : edges[cand_out]) {
                    if (!in_nei_c[x]) ++new_size_outer;
                }
                for (auto x : edges[cand_in]) {
                    ++in_nei_c[x];
                }
            }
            double p = (nodes[cand_out] * cur_size_outer) / (nodes[cand_in] * new_size_outer);
            if (gen_p >= p)
                return false;
            inner_update(cand_out, false);
            if (!is_connected(cand_in)) {
                inner_update(cand_out, true);
                return false;
            }
            inner_update(cand_in, true);
            update_outer_nodes(cand_in, cand_out);
            return true;
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
                if (!is_connected(cand)) {
                    return false;
                }
                inner_update(cand, erase);
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
