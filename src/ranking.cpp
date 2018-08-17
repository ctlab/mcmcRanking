#include <Rcpp.h>
#include <vector>
#include "utils.h"
#include "mcmc.h"
#include "hsa.h"

using namespace Rcpp;
using namespace std;
using namespace mcmc;


void
dfs(unsigned i, vector <vector<unsigned>> &edges, vector<unsigned> &new_comp, vector<char> &used, HSA &not_ranked) {
    used[not_ranked.get_index(i)] = true;
    new_comp.push_back(i);
    for (unsigned x : edges[i]) {
        size_t ind_x = not_ranked.get_index(x);
        if (ind_x == -1)
            continue;
        if (!used[ind_x])
            dfs(x, edges, new_comp, used, not_ranked);
    }
}


void
cut_points_dfs(int v, int p, int &timer, vector <vector<unsigned>> &edges, vector<unsigned> &ranked, vector<char> &used,
               int tin[], int fup[], vector<char> &is_cut_point) {
    used[v] = true;
    tin[v] = timer++;
    fup[v] = timer;
    int children = 0;
    for (unsigned to : edges[v]) {
        if (to == p || ranked[to])
            continue;
        if (used[to]) {
            fup[v] = min(fup[v], tin[to]);
        } else {
            cut_points_dfs(to, v, timer, edges, ranked, used, tin, fup, is_cut_point);
            fup[v] = min(fup[v], fup[to]);
            if (fup[to] >= tin[v] && p != -1)
                is_cut_point[v] = true;
            ++children;
        }
    }
    if (p == -1 && children > 1)
        is_cut_point[v] = true;
}

vector<char> cut_points_dfs(int v, vector <vector<unsigned>> &edges, vector<unsigned> &ranked) {
    size_t n = edges.size();
    vector<char> is_cut_point(n, false);
    vector<char> used(n, false);
    int tin[n];
    int fup[n];
    int timer = 0;
    cut_points_dfs(v, -1, timer, edges, ranked, used, tin, fup, is_cut_point);
    return is_cut_point;
}

bool
get_best_comp(vector <vector<unsigned>> &edges, vector<double> &nodes, HSA &not_ranked,
              vector<unsigned> &comp, vector<unsigned> &to_remove, int i) {
    vector<char> used(not_ranked.size(), false);
    used[not_ranked.get_index(i)] = true;
    to_remove.push_back(i);
    for (int j = 0; j < nodes.size(); ++j) {
        size_t ind_j = not_ranked.get_index(j);
        if (ind_j == -1)
            continue;
        if (!used[ind_j]) {
            vector<unsigned> new_comp;
            dfs(j, edges, new_comp, used, not_ranked);
            if (new_comp.size() > comp.size())
                comp.swap(new_comp);
            for (unsigned x : new_comp) {
                to_remove.push_back(x);
                if (nodes[x] > nodes[i])
                    return true;
            }
        }
    }
    return false;
}

// [[Rcpp::export]]
IntegerVector probabilistic_rank_internal(IntegerMatrix edgelist, DataFrame df_nodes) {
    NumericVector q = df_nodes["q"];
    unsigned n = q.size();
    vector<double> nodes(q.begin(), q.end());
    vector <vector<unsigned>> edges = adj_list(edgelist, n);

    vector<unsigned> ranked(n, 0);

    HSA not_ranked(n);
    for (size_t i = 0; i < n; ++i) {
        not_ranked.insert(i);
    }

    double sum_p = 0;
    for (double x : nodes)
        sum_p += 1 - x;
    double sum_q = n - sum_p;

    unsigned cand_remove = 0;
    while (not_ranked.size() != 0) {
        double best_q = sum_q;
        double best_p = -1;

        size_t nr_v = not_ranked.get(0);
        vector<char> is_cut_point = cut_points_dfs(nr_v, edges, ranked);

        size_t ind = nr_v;
        for (size_t i : not_ranked.get_all()) {
            if (nodes[i] >= nodes[ind]) {
                if (nodes[i] == nodes[ind] && is_cut_point[i]) {
                    continue;
                }
                ind = i;
            }
        }
        if (!is_cut_point[ind]) {
            ranked[ind] = not_ranked.size();
            not_ranked.erase(ind);
            sum_p -= 1 - nodes[ind];
            sum_q -= nodes[ind];
            continue;
        }
        for (size_t i : not_ranked.get_all()) {
            double cur_p;
            double cur_q;
            if (is_cut_point[i]) {
                vector<unsigned> comp;
                vector<unsigned> to_remove;
                if (get_best_comp(edges, nodes, not_ranked, comp, to_remove, i))
                    continue;
                cur_p = 0;
                for (unsigned x : comp)
                    cur_p += 1 - nodes[x];
                cur_q = comp.size() - cur_p;
            } else {
                cur_p = sum_p - (1 - nodes[i]);
                cur_q = sum_q - nodes[i];
            }
            double sin = (best_p - sum_p) * (cur_q - sum_q) - (best_q - sum_q) * (cur_p - sum_p);
            bool smaller = cur_p * cur_p + cur_q * cur_q < best_p * best_p + best_q * best_q;
            if (sin > 0 || sin == 0 && smaller || best_p == -1) {
                best_p = cur_p;
                best_q = cur_q;
                cand_remove = i;
            }
        }
        if (is_cut_point[cand_remove]) {
            vector<unsigned> comp;
            vector<unsigned> to_remove;
            get_best_comp(edges, nodes, not_ranked, comp, to_remove, cand_remove);
            unsigned rank_all = not_ranked.size();
            for (unsigned x : to_remove) {
                ranked[x] = rank_all;
                not_ranked.erase(x);
            }
        } else {
            ranked[cand_remove] = not_ranked.size();
            not_ranked.erase(cand_remove);
        }
        sum_p = best_p;
        sum_q = best_q;
    }
    return IntegerVector(ranked.begin(), ranked.end());
}
