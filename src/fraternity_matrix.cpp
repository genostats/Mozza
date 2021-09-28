#include <Rcpp.h>
#include "mozza.h"
using namespace Rcpp;

namespace mozza {
NumericMatrix fraternity_matrix(std::vector<zygote> & ZYG) {
  int n = ZYG.size();
  double total_length(0);
  if(n > 0)
    total_length = ZYG[0].first.genome_length; // on suppose que tout a la même longueur partout...
  NumericMatrix K(n,n);
  for(int i = 0; i < n; i++) {
    for(int j = i; j < n; j++) {
      auto IBD = IBD_length(ZYG[i], ZYG[j]);
      K(j,i) = (std::get<2>(IBD)) / total_length;
    }
  }
  // symetriser
  for(int i = 1; i < n; i++) 
    for(int j = 0; j < i; j++)
      K(j,i) = K(i,j);
  return K;
}

} // end namespace
