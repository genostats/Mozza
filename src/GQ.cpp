#include <Rcpp.h>
using namespace Rcpp;

// remplace GQ <- round( apply(phtab, 1, median) ); GQ[GQ > 99] <- 99
// phtab etant une matrice trois colonnes
// [[Rcpp::export]]
IntegerVector GQ(NumericMatrix phtab) {
  int n = phtab.nrow();
  IntegerVector R = no_init(n);
  for(int i = 0; i < n; i++) {
    int a0 = (int) round(phtab(i,0)); 
    int a1 = (int) round(phtab(i,1));
    int a2 = (int) round(phtab(i,2));
    a0 = a0>99?99:a0;
    a1 = a1>99?99:a1;
    a2 = a2>99?99:a2;
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

