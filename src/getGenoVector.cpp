#include <Rcpp.h>
#include "gaston/matrix4.h"

// [[Rcpp::export]]
IntegerVector getGenoVector(XPtr<matrix4> pA, unsigned int snpIndex, bool Rindex = true) {
  if(Rindex) snpIndex--;
  if((snpIndex < 0) | (snpIndex >= pA->nrow)) stop("bad index");
  IntegerVector X(pA->ncol);
  int tab[4] = {0, 1, 2, NA_INTEGER};
  for(size_t i = 0; i < pA->true_ncol-1; i++) { 
    size_t x = pA->data[snpIndex][i];
    for(unsigned int ss = 0; ss < 4; ss++) {
      X[4*i+ss] = tab[ x&3 ]; // ((x&3) != 3)?(x&3):NA_INTEGER;
      x >>= 2;
    }
  }
  size_t i = pA->true_ncol-1;
  size_t x = pA->data[snpIndex][i];
  for(unsigned int ss = 0; ss < 4 && 4*i+ss < pA->ncol; ss++) {
    X[4*i+ss] = tab[ x&3 ]; // ((x&3) != 3)?(x&3):NA_INTEGER;
    x >>= 2;
  }
  return X;
}

