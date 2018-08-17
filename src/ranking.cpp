#include <Rcpp.h>
#include <vector>
#include "utils.h"
#include "mcmc.h"

using namespace Rcpp;
using namespace std;
using mcmc::Graph;


void dfs(unsigned i, vector<vector<unsigned>> &edges, vector<unsigned> &new_comp, vector<bool> &used) {
    used[i] = true;
    new_comp.push_back(i);
    for (unsigned x : edges[i])
        if (!used[x])
            dfs(x, edges, new_comp, used);
}


void cut_points_dfs(int v, int p, int &timer, vector<vector<unsigned>> &edges, vector<unsigned> &ranked, bool used[],
                    int tin[], int fup[], bool is_cut_point[]) {
    used[v] = true;
    tin[v] = timer++;
    fup[v] = timer;
    int children = 0;
    for (size_t i = 0; i < edges[v].size(); ++i) {
        int to = edges[v][i];
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


bool
get_best_comp(vector<vector<unsigned>> &edges, vector<double> &nodes, vector<unsigned> &ranked, vector<unsigned> &comp,
              vector<unsigned> &to_remove, int i) {
    vector<bool> used(ranked.begin(), ranked.end());
    used[i] = true;
    to_remove.push_back(i);
    for (int j = 0; j < nodes.size(); ++j) {
        if (!used[j]) {
            vector<unsigned> new_comp;
            dfs(j, edges, new_comp, used);
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
    vector<vector<unsigned>> edges = adj_list(edgelist, n);

    vector<unsigned> ranked(n, 0);
    unsigned not_ranked_n = n;

    double sum_p = 0;
    for (double x : nodes)
        sum_p += 1 - x;
    double sum_q = n - sum_p;

    unsigned cand_remove = 0;
    while (not_ranked_n != 0) {
        double best_q = sum_q;
        double best_p = -1;

        bool used[n];
        bool is_cut_point[n];
        memset(used, 0, sizeof used);
        memset(is_cut_point, 0, sizeof is_cut_point);
        int tin[n];
        int fup[n];
        int timer = 0;
        for (int i : edges[cand_remove]) {
            if (!ranked[i]) {
                cut_points_dfs(i, -1, timer, edges, ranked, used, tin, fup, is_cut_point);
                break;
            }
        }
        size_t ind = -1;
        for(int i = 0; i < n; ++i){
            if(ranked[i])
                continue;
            if(ind == -1) {
                ind = i;
            } else if(nodes[i] >= nodes[ind]) {
                if (nodes[i] == nodes[ind] && is_cut_point[i]) {
                    continue;
                }
                ind = i;
            }
        }
        if(!is_cut_point[ind]){
            ranked[ind] = not_ranked_n--;
            sum_p = sum_p - (1 - nodes[ind]);
            sum_q =  sum_q - nodes[ind];
            continue;
        }
        for (int i = 0; i < n; ++i) {
            if (ranked[i])
                continue;
            double cur_p;
            double cur_q;
            if (is_cut_point[i]) {
                vector<unsigned> comp;
                vector<unsigned> to_remove;
                if (get_best_comp(edges, nodes, ranked, comp, to_remove, i))
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
            get_best_comp(edges, nodes, ranked, comp, to_remove, cand_remove);
            for (unsigned x : to_remove) {
                ranked[x] = not_ranked_n;
            }
            not_ranked_n -= to_remove.size();
        } else {
            ranked[cand_remove] = not_ranked_n--;
        }
        sum_p = best_p;
        sum_q = best_q;
    }
    return IntegerVector(ranked.begin(), ranked.end());
}
