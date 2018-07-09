#include <Rcpp.h>
#include <vector>
#include <queue>
#include "utils.h"
using namespace Rcpp;
using namespace std;

bool is_connected(vector<vector<unsigned>> edges, bool from_inner[]) {
  int n = edges.size();
  vector<bool> used(n, false);
  queue<unsigned> q;
  for(int i = 0; i < n; ++i){
    if(from_inner[i]){
      used[i] = true;
      q.push(i);
      break;
    }else if(i == n - 1){
      return true;
    }
  }
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
  for(int i = 0; i < n; ++i){
    if(from_inner[i] && !used[i]){
      return false;
    }
  }
  return true;
}

// [[Rcpp::export]]
NumericVector real_prob_internal(IntegerMatrix edgelist, NumericVector likelihood){
  vector<vector<unsigned>> edges = adj_list(edgelist, likelihood.size());
  vector<double> scores(likelihood.size(), 0.0);
  double sumscores = 0;
  int n = likelihood.size();
  bool x[n];
  for(int i = 1; i < n + 1; ++i){
    for(int j = 0; j < n; ++j)  x[j] = i + j < n ? false : true;
    do{
      if(is_connected(edges, x)){
        double score = 1;
        for(int i = 0; i < n; ++i) if(x[i]) score *= likelihood[i];
        sumscores += score;
        for(int i = 0; i < n; ++i) if(x[i]) scores[i] += score;
      }
    } while (next_permutation(x, x + n));
  }
  for(int i = 0; i < n; ++i)
    scores[i] /= sumscores;
  return NumericVector (scores.begin(), scores.end());
}

//' Accurate sum of numbers
//'
//' Sums every time only two smallest numbers.
//' This allows calculate the sum of numbers accurate enough.
//' @export
// [[Rcpp::export]]
NumericVector accurate_sum(NumericVector v){
  int n = v.size();
  if(n == 0){
    return NumericVector::create(0);
  }
  if(n == 1){
    return NumericVector::create(v[0]);
  }
  vector<double> a(n);
  copy(v.begin(), v.end(), a.begin());
  sort(a.begin(), a.end());
  a[0] = a[0] + a[1];
  int i = 2;
  int j = 0;
  int end = 1;
  while (i < n || j + 1 != end){
    if(i < n){
      if(a[i] <= a[j]){
        if(i + 1 < n){
          if(a[i+1] <= a[j]){
            a[end++] = a[i++] + a[i++];
            continue;
          }
        }
        a[end++] = a[i++] + a[j++];
      }else if(j + 1 < end && a[j+1] < a[i]){
        a[end++] = a[j++] + a[j++];
      }else{
        a[end++] = a[i++] + a[j++];
      }
    }else{
      a[end++] = a[j++] + a[j++];
    }
  }
  return NumericVector::create(a[j]);
}
