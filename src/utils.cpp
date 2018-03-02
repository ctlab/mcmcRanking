#include <Rcpp.h>
#include <vector>
using namespace Rcpp;
using namespace std;

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
