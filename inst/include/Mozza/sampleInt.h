#include <Rcpp.h>
#ifndef __sampleInt__
#define __sampleInt__
template<typename DV>
inline int sampleInt(DV & probs) {
  double u = R::unif_rand();
  int n = probs.size();
  for(int i = 0; i < n; i++) {
    auto p = probs[i];
    if(u < p) return i;
    u -= p;
  }
  Rcpp::stop("In sampleInt, probs doesn't sum to 1");
}
#endif

