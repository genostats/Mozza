#include <Rcpp.h>
#include "mozza.h"
using namespace Rcpp;

namespace mozza {
NumericMatrix kinship_matrix(std::vector<zygote> & ZYG) {
  int n = ZYG.size();
  double total_length(0);
  if(n > 0)
    total_length = ZYG[0].first.genome_length; // on suppose que tout a la même longueur partout...
  NumericMatrix K(n,n);
  for(int i = 0; i < n; i++) {
    // le coeff diagonal
    double HBD = HBD_length(ZYG[i]);
    K(i,i) = 1.0 + HBD / total_length; // 1 + coeff de consanguinité
    for(int j = i+1; j < n; j++) {
      // auto IBD = IBD_length(ZYG[i], ZYG[j]);
      // K(j,i) = (0.5*std::get<1>(IBD) + std::get<2>(IBD)) / total_length;
      double R = relatednessLength(ZYG[i], ZYG[j]);
      K(j,i) = 2.0 * R / total_length;
    }
  }
  // symetriser
  for(int i = 1; i < n; i++) 
    for(int j = 0; j < i; j++)
      K(j,i) = K(i,j);
  return K;
}

} // end namespace
