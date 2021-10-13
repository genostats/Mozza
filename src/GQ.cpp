#include <Rcpp.h>
using namespace Rcpp;

// remplace apply(phtab, 1, median)
// phtab etant une matrice trois colonnes
// [[Rcpp::export]]
NumericVector GQ(NumericMatrix phtab) {
  int n = phtab.nrow();
  NumericVector R = no_init(n);
  for(int i = 0; i < n; i++) {
    double a0 = phtab(i,0); 
    double a1 = phtab(i,1);
    double a2 = phtab(i,2);
    // c'est moche mais trois tests maxi...
    if(a0 <= a1) {
      if(a1 <= a2) { // a0 <= a1 <= a2
        R[i] = a1; continue;
      }
      if(a2 <= a0) { // a2 <= a0 <= a1
        R[i] = a0; continue;
      } 
      R[i] = a2; continue; // a0 < a2 < a1
    }
    if(a1 >= a2) {
      R[i] = a1; continue; // a2 <= a1 < a0
    }
    if(a0 <= a2) {
      R[i] = a0; continue; // a1 < a0 <= a2
    }
    R[i] = a2; // a1 < a2 < a0
  }
  return R;
}

