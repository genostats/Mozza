#include <Rcpp.h>
using namespace Rcpp;


// remplace apply (phtab, 1,  function (x) which (x == 0)[1]))
// phtab etant une matrice trois colonnes
// [[Rcpp::export]]
IntegerVector GT(NumericMatrix phtab) {
  int n = phtab.nrow();
  IntegerVector R = no_init(n);
  for(int i = 0; i < n; i++) {
    double a0 = phtab(i,0); 
    double a1 = phtab(i,1);
    double a2 = phtab(i,2);
    if(a0 == 0) { 
      R[i] = 1; continue;
    }
    if(a1 == 0) {
      R[i] = 2; continue;
    }
    if(a2 == 0) {
      R[i] = 3; continue;
    }
    stop("Bad phred matrix");
  }
  return R;
}
