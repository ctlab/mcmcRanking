#ifndef UTILS_MCMC_H
#define UTILS_MCMC_H

#include <Rcpp.h>
#include <vector>

using namespace Rcpp;
using namespace std;

vector<vector<unsigned>> adj_list(IntegerMatrix edgelist, size_t gorder);

#endif
