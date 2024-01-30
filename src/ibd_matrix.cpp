#include <Rcpp.h>
#include "Mozza.h"
using namespace Rcpp;

namespace mozza {
NumericMatrix ibd_matrix(std::vector<mosaic> & HAP) {
  int n = HAP.size();
  double total_length(0);
  if(n > 0)
    total_length = HAP[0].genome_length; // on suppose que tout a la même longueur partout...
  NumericMatrix K(n,n);
  for(int i = 0; i < n; i++) {
    // le coeff diagonal
    K(i,i) = 1.0;
    for(int j = i+1; j < n; j++) {
      double R = IBD_sharing(HAP[i], HAP[j]);
      K(j,i) = R / total_length;
    }
  }
  // symetriser
  for(int i = 1; i < n; i++) 
    for(int j = 0; j < i; j++)
      K(j,i) = K(i,j);
  return K;
}

} // end namespace
