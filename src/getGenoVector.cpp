#include <Rcpp.h>
#include "gaston/matrix4.h"

// [[Rcpp::export]]
IntegerVector getGenoVector(XPtr<matrix4> pA, int snpIndex, bool Rindex = true) {
  if(Rindex) snpIndex--;
  if((snpIndex < 0) | (snpIndex >= pA->nrow)) stop("bad index");
  IntegerVector X(pA->ncol);
  int tab[4] = {0, 1, 2, NA_INTEGER};
  for(int i = 0; i < pA->true_ncol-1; i++) { 
    uint8_t x = pA->data[snpIndex][i];
    for(int ss = 0; ss < 4; ss++) {
      X[4*i+ss] = tab[ x&3 ]; // ((x&3) != 3)?(x&3):NA_INTEGER;
      x >>= 2;
    }
  }
  int i = pA->true_ncol-1;
  uint8_t x = pA->data[snpIndex][i];
  for(int ss = 0; ss < 4 && 4*i+ss < pA->ncol; ss++) {
    X[4*i+ss] = tab[ x&3 ]; // ((x&3) != 3)?(x&3):NA_INTEGER;
    x >>= 2;
  }
  return X;
}

